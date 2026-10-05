# VariantRail

VariantRail is a small backend for a genome variant annotation pipeline with a
reproducible provenance trail. It ingests a VCF subset, filters records by
quality and depth, annotates every ALT allele with a built-in gene/consequence
table, summarises the result, and records one SHA-256 link per pipeline step so
any run can be re-derived and checked byte for byte.

The initial release intentionally supports a compact public contract:

- a sample is VCF text that is parsed and validated at upload time;
- a run replays four fixed steps: `ingest`, `filter`, `annotate`, `summarize`;
- each step stores the hash of its input, the hash of its output, and the
  previous link's hash, so the trail is a chain;
- hashes depend only on content and parameters, never on the sample id, the
  run id, wall-clock time, or process state;
- the annotation table is built in and its content hash is itself a
  provenance input;
- every output variant keeps the original VCF line verbatim.

## Requirements

- Python 3.11 or newer
- no third-party runtime dependencies

## Run the service

```bash
PYTHONPATH=src python -m variantrail.server --host 127.0.0.1 --port 18083 --database variantrail.db
```

The process prints `VariantRail listening on http://127.0.0.1:18083` after it
has bound the port.

## HTTP API

All request and response bodies are JSON except the `jsonl`, `tsv` and `vcf`
run export endpoints, and unknown fields are rejected. Both `POST` endpoints
require an `Idempotency-Key` header; repeating a key returns the stored first
response verbatim, and reusing a key for another operation is a `conflict`.

### Health

`GET /health` returns `{"status":"ok"}`.

### Create a sample

```http
POST /samples
Idempotency-Key: sample-1
Content-Type: application/json

{"id":"trio-1","vcf":"##fileformat=VCFv4.2\n##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t11856378\trs1\tG\tA\t60\tPASS\tDP=42\nchr1\t11856378\trs2\tG\tT\t20\tPASS\tDP=42\nchr7\t55019017\t.\tG\tGA\t50\tPASS\tDP=31\nchr17\t43093445\t.\tC\tT\t99\tPASS\tDP=55\n"}
```

Returns HTTP 201 with a deterministic projection of the parsed sample:

```json
{"chromosomes":["chr1","chr7","chr17"],"id":"trio-1","meta_lines":2,"records":4,
 "sha256":"45520c9b2b4798adafe3bdbaf137c567536527581c259dc2d92c391baf2a94ef"}
```

`sha256` is the SHA-256 of the uploaded VCF text. `id` must match
`[A-Za-z0-9][A-Za-z0-9._:-]{0,99}` and is unique (`conflict` otherwise).

### The supported VCF subset

The text is split on `\n`; one trailing newline is optional and a single
trailing `\r` per line is ignored, so CRLF files are accepted. Every violation
raises `validation_error` with a message shaped `line <n>: ...`.

1. The first line must be `##fileformat=VCFv4.<digits>`.
2. Zero or more `##key=value` meta lines follow, `key` matching
   `[A-Za-z][A-Za-z0-9_.-]*` and `value` non-empty.
3. Exactly one `#CHROM` line follows whose tab-separated columns are exactly
   `#CHROM POS ID REF ALT QUAL FILTER INFO`. `FORMAT` and sample columns are
   not supported and are rejected.
4. One or more data lines follow, each with exactly 8 tab-separated fields:
   - `CHROM` matches `[A-Za-z0-9_.]+`; chromosome blocks must be contiguous, so
     a chromosome may not reappear once another chromosome started;
   - `POS` is a positive integer and non-decreasing within a chromosome, and
     `(CHROM, POS, REF, ALT)` is unique in the file;
   - `ID` is `.` (becomes `null`) or a value without whitespace, `;` or `,`;
   - `REF` is a non-empty `[ACGTN]+` sequence;
   - `ALT` is one or more `[ACGTN]+` sequences separated by `,` without
     repeats; symbolic alleles such as `<DEL>` or `*` are rejected;
   - `REF` and `ALT` are upper-cased before storage;
   - `QUAL` is `.` (becomes `null`) or a non-negative decimal number;
   - `FILTER` is `.` (empty list), `PASS`, or `;`-separated tokens matching
     `[A-Za-z0-9_.]+`; `PASS` must not be combined with another token;
   - `INFO` is `.` (empty object) or `;`-separated `KEY` / `KEY=VALUE` entries
     with `KEY` matching `[A-Za-z_][A-Za-z0-9_.]*`, no repeated keys and a
     non-empty whitespace-free value; `DP`, when present, must be a
     non-negative integer.

Blank lines, header lines after the first data record, a missing `#CHROM`
header, and a file with no data records are rejected.

### Start a run

```http
POST /samples/trio-1/runs
Idempotency-Key: run-1
Content-Type: application/json

{"id":"run-1","params":{"min_qual":30,"min_dp":10,"pass_only":true,"genes":["MTHFR","EGFR","BRCA1"]}}
```

`params` is optional; the body must contain `id` and may contain only `params`
besides it. The run executes synchronously, so a returned run is always
`succeeded`. Returns HTTP 201:

```json
{"id":"run-1","sample_id":"trio-1",
 "sample_sha256":"45520c9b2b4798adafe3bdbaf137c567536527581c259dc2d92c391baf2a94ef",
 "params":{"genes":["BRCA1","EGFR","MTHFR"],"impacts":[],"min_dp":10,"min_qual":30.0,"pass_only":true},
 "status":"succeeded","provenance_length":4,
 "provenance_head":"40807a2cebf5110a5a3be231ab78ca690823b9def8414dd373314354cac65b40",
 "statistics":{"records_total":4,"records_kept":3,"records_filtered":1,
   "filter_counts":{"allele_filtered":0,"dp_below_min":0,"not_pass":0,"qual_below_min":1},
   "alleles_total":4,"alleles_kept":3,
   "consequence_counts":{"frameshift_variant":1,"missense_variant":1,"stop_gained":1},
   "impact_counts":{"HIGH":2,"MODERATE":1},
   "gene_counts":{"BRCA1":1,"EGFR":1,"MTHFR":1}}}
```

### Run parameters

Every key is optional; the stored `params` object always contains all five
core keys, plus `record_filter` only when one was supplied. `genes` and
`impacts` are stored sorted, so request order never changes a hash.

- `pass_only` (boolean, default `false`) — keep only records whose `FILTER` is
  exactly `["PASS"]`.
- `min_qual` (non-negative finite number, default `0.0`, stored as a float) —
  keep records with `QUAL >= min_qual`.
- `min_dp` (non-negative integer, default `0`) — keep records with
  `DP >= min_dp`.
- `genes` (array of non-empty strings, default `[]` = no restriction) — after
  annotation keep only alleles whose gene is listed. An allele with no table
  match has gene `null` and is dropped when `genes` is non-empty.
- `impacts` (array, default `[]` = no restriction) — after annotation keep only
  alleles whose impact is listed; values must come from `HIGH`, `MODERATE`,
  `LOW`, `MODIFIER`, `UNKNOWN`.
- `record_filter` (expression object, optional) — a structured, nested
  record-level predicate applied after `pass_only`, `min_qual` and `min_dp` and
  before the gene/impact allowlists. When omitted the stored `params` object
  holds only the other five keys and the run is byte-for-byte identical to an
  older run; when present it is stored, normalized, under `record_filter`.

A record that passed the record-level filters but whose alleles are all removed
by `genes`/`impacts` is dropped as well.

### Record filter expressions

Every expression node is either a logical node or a condition object; mixing
both forms in one node is rejected:

- `{"all":[<expr>, ...]}` — true when every child is true; the array is
  non-empty.
- `{"any":[<expr>, ...]}` — true when at least one child is true; the array is
  non-empty.
- `{"not":<expr>}` — true when the single child expression is false.
- `{"field":"<field>","op":"<op>","value":<value>}` — a condition; `value` is
  required for every op except `exists`, which must omit it.

`field` is one of `chrom`, `pos`, `id`, `ref`, `alt`, `qual`, `dp`, `filter`
or `info.<KEY>` (any INFO key). `op` is one of:

- `eq` — direct scalar comparison. For `alt` and `filter` the condition is
  true when any element equals the value.
- `in` — `value` must be a non-empty array; true on membership. For `alt` and
  `filter` it is true when any record element appears in the array.
- `lt`, `lte`, `gt`, `gte` — numeric ordering; `value` must be a number and
  the field must be `pos`, `qual`, `dp`, or an `info.<KEY>` whose stored value
  parses fully to a finite decimal number.
- `exists` — takes no `value`; true when the field is present and non-null. An
  INFO bare flag counts as present.

A missing field, a `null` value (`ID`/`QUAL` `.`), or an INFO value that does
not parse as a finite decimal makes every op other than `exists` false; wrap
such a condition in `not` to test for absence. The stored expression is
normalized (`in` arrays are sorted, numeric values stored as floats), so
equivalent submissions hash identically. Expressions may nest at most 16
levels. Any unknown field or op, malformed node, empty `all`/`any`,
field/op type mismatch, `exists` with a `value`, another op without one, or
nesting beyond 16 levels makes run creation return 400 `validation_error`
without saving the run or the idempotency response.

### Filter accounting and statistics

`filter_counts` attributes every filtered record to exactly one reason, in this
precedence: `not_pass` (`pass_only` set and `FILTER` is not `["PASS"]`);
`qual_below_min` (`QUAL` below `min_qual`, or `null` while `min_qual > 0`);
`dp_below_min` (`DP` below `min_dp`, or `null` while `min_dp > 0`);
`rule_mismatch` (a `record_filter` was supplied and its expression evaluated
false for the record — present in the run statistics and annotate-step counts,
even when zero); `allele_filtered` (survived the record-level filters but no
allele survived the allowlists). Hence
`records_filtered == sum(filter_counts.values())` and
`records_kept == records_total - records_filtered` always hold.

`alleles_total` counts every ALT allele of every record and `alleles_kept`
counts the alleles present in the run's variants. `gene_counts` uses the key
`NA` for alleles with no table match. `consequence_counts`, `impact_counts` and
`gene_counts` omit zero counts and are sorted by key.

### Get a run

`GET /runs/run-1` returns the run document shown above.

### Get run variants

`GET /runs/run-1/variants` returns `{"count":<n>,"run_id":"run-1","variants":[...]}`
in file order. Each variant repeats the parsed fields plus `raw`, the original
line, so every output stays traceable to the input:

```json
{"alts":["A"],
 "annotations":[{"allele":"A","allele_index":1,"consequence":"missense_variant","gene":"MTHFR","impact":"MODERATE"}],
 "chrom":"chr1","dp":42,"filter":["PASS"],"id":"rs1","info":{"DP":"42"},"line":4,
 "pos":11856378,"qual":60.0,"raw":"chr1\t11856378\trs1\tG\tA\t60\tPASS\tDP=42","ref":"G"}
```

`line` is the 1-based line of the record in the uploaded VCF and `allele_index`
is the 1-based VCF ALT index.

### Export run results

`GET /runs/run-1/exports/jsonl`, `GET /runs/run-1/exports/tsv` and
`GET /runs/run-1/exports/vcf` hand the
retained variants of a saved run to downstream analysis as a stable payload.
The body is generated on demand from the stored run; it is never written to
disk, and the run, its variants and its provenance are not modified.

`jsonl` is served as `application/x-ndjson; charset=utf-8`: one variant object
per line in the original file order, using the same fields and values as
`GET /runs/{run_id}/variants`. Each line is UTF-8 JSON with no indentation and
keys sorted lexicographically, terminated by a single LF; the body ends with
the last line's LF (no trailing blank line). A run with no variants returns an
empty body.

`tsv` is served as `text/tab-separated-values; charset=utf-8` and always opens
with the fixed header:

```
chrom	pos	id	ref	alt	allele_index	gene	consequence	impact	qual	dp	filter	info	line
```

Each retained annotation is one row, so a record with several retained ALT
alleles spans several rows; variant file order and annotation order are kept.
Formatting rules:

- `id`, `gene`, `qual` or `dp` is `NA` when null; `qual` otherwise keeps the
  JSON decimal form (e.g. `60.0`);
- an empty `FILTER` is `.`; several filters are joined with `;`;
- `INFO` is `;`-joined `KEY=VALUE` pairs sorted by key; bare flags are written
  as `KEY`; empty INFO is `.`;
- numbers use the JSON decimal representation and text is not quoted.

A run with no variants returns just the header line. An unknown `run_id`
returns 404 `not_found`; any `format` other than `jsonl`, `tsv` or `vcf`
returns 400 `validation_error`, both with the standard error object.

`vcf` is served as `text/x-variant-call-format; charset=utf-8`: the retained
variants rendered back as VCF v4.2 text (UTF-8, LF line endings). The header
opens with `##fileformat=VCFv4.2`, repeats the uploaded sample's other meta
lines in their original order, adds `##source=variantrail`, and declares four
export INFO keys — `ALLELE_INDEX` (Integer), `VT_GENE`, `VT_CONSEQUENCE` and
`VT_IMPACT` (String, all `Number=1`) — before the tab-separated
`#CHROM POS ID REF ALT QUAL FILTER INFO` column header. Each retained ALT
allele is one record, in variant file order and `allele_index` order, with
`ALT` holding only the current allele. `ID`, `QUAL` and an empty `FILTER` are
`.` when missing; `INFO` lists the original entries sorted by key, then the
four export keys (`VT_GENE` is `NA` when no gene is annotated). If an original
INFO field uses `ALLELE_INDEX` or one of the `VT_*` keys, the export returns
400 `validation_error` naming the conflicting key. A run with no variants
returns just the header, and repeated requests are byte-identical.

### Compare two runs

`GET /runs/{run_id}/compare/{other_run_id}` compares the retained ALT alleles
of two finished runs — `run_id` on the left, `other_run_id` on the right. It
only reads the stored runs: no pipeline step is re-executed and neither run,
its variants, its exports, nor its provenance is modified. The runs may belong
to different samples. No `Idempotency-Key` is required.

An allele's identity is the tuple `(chrom, pos, ref, alt)`; duplicates of the
same identity within one run count once. Returns HTTP 200:

```json
{"left":{"run_id":"run-1","sample_id":"trio-1","sample_sha256":"...","allele_count":8},
 "right":{"run_id":"run-2","sample_id":"trio-1","sample_sha256":"...","allele_count":4},
 "counts":{"shared":4,"left_only":4,"right_only":0,"union":8},
 "shared":[{"chrom":"chr1","pos":11856378,"ref":"G","alt":"T",
   "gene":"MTHFR","consequence":"stop_gained","impact":"HIGH"}],
 "left_only":[...],
 "right_only":[...],
 "gene_summary":{"MTHFR":{"shared":1,"left_only":1,"right_only":0}},
 "impact_summary":{"HIGH":{"shared":4,"left_only":0,"right_only":0}}}
```

- `shared`, `left_only` and `right_only` hold every distinct allele on both
  sides, on the left only, and on the right only; each element carries
  `chrom`, `pos`, `ref`, `alt`, `gene`, `consequence` and `impact`, with
  `gene` kept as JSON `null` when the allele has no table match. Comparing a
  run with itself puts every allele in `shared`.
- Each array is sorted by `chrom` (text order), then `pos` (numeric), then
  `ref`, then `alt`, so repeated requests are byte-for-byte identical.
- `counts` mirrors the three array lengths plus `union`, their total.
- `gene_summary` and `impact_summary` count the alleles of each bucket per
  gene and per impact; keys are sorted lexicographically, a `null` gene
  appears as `NA`, and impacts keep their original values.

An unknown `run_id` or `other_run_id` returns 404 `not_found` with the
message `run <run_id> was not found`.

### Compare a cohort of runs

`GET /runs/{run_id}/cohort/{other_run_ids}` compares the retained ALT alleles
of two or more finished runs at once. `other_run_ids` is two or more run ids
joined with commas; the path order (the path `run_id` first, then the comma
list) is the cohort queue order. Like the two-run comparison it only reads the
stored runs: no pipeline step is re-executed and no run, variant, export or
provenance record is modified. No `Idempotency-Key` is required.

An allele's identity is again `(chrom, pos, ref, alt)`; duplicates of the same
identity within one run count once, and annotations are taken as stored — no
filtering or re-annotation happens. Returns HTTP 200:

```json
{"run_count":3,
 "runs":[{"run_id":"run-1","sample_id":"trio-1","sample_sha256":"...","allele_count":8},
         {"run_id":"run-2","sample_id":"trio-1","sample_sha256":"...","allele_count":4},
         {"run_id":"run-3","sample_id":"trio-2","sample_sha256":"...","allele_count":3}],
 "alleles":[{"chrom":"chr1","pos":11856378,"ref":"G","alt":"A","gene":"MTHFR",
   "consequence":"missense_variant","impact":"MODERATE","run_count":2,
   "present_in":["run-1","run-3"]}],
 "counts":{"unique_alleles":9,"core_alleles":2,"variable_alleles":3,
   "private_alleles":{"run-1":3,"run-2":1,"run-3":0}},
 "frequency_summary":{"1":4,"2":3,"3":2},
 "gene_summary":{"BRCA1":{"unique_alleles":1,"core_alleles":1,"variable_alleles":0},
                "MTHFR":{"unique_alleles":2,"core_alleles":0,"variable_alleles":2}},
 "impact_summary":{"HIGH":{"unique_alleles":4,"core_alleles":2,"variable_alleles":2}}}
```

- `runs` repeats the queued runs in queue order; each entry has `run_id`,
  `sample_id`, `sample_sha256` and `allele_count`, the number of distinct
  retained alleles in that run.
- `alleles` holds every distinct allele of the union (different ALTs at the
  same coordinate are separate entries). Each element carries `chrom`, `pos`,
  `ref`, `alt`, `gene` (`null` when the allele has no table match),
  `consequence`, `impact`, `run_count` — the number of cohort runs containing
  it — and `present_in`, the run ids of those runs in queue order. The array
  is sorted by `chrom` (text order), then `pos` (numeric), then `ref`, then
  `alt`, so repeated requests are byte-for-byte identical.
- `counts.unique_alleles` is the total number of entries in `alleles`;
  `core_alleles` counts alleles present in all `run_count` runs,
  `variable_alleles` alleles present in between 2 and `run_count - 1` runs,
  and `private_alleles` breaks the alleles present in exactly one run down by
  `run_id`, in queue order (runs without a private allele show `0`).
- `frequency_summary` has the string keys `"1"` through
  `"<run_count>"` and counts alleles present in exactly that many runs, with
  missing frequencies filled in as `0`.
- `gene_summary` and `impact_summary` group `unique_alleles`, `core_alleles`
  and `variable_alleles` per gene and per impact; a `null` gene appears as
  `NA`, categories with no alleles are omitted, and keys are sorted
  lexicographically. Private alleles are not separately listed there
  (`unique = core + variable + private`).

Changing the queue order returns the same data apart from the order of `runs`
and of each `present_in`; the allele set, counts and summaries are unchanged,
and repeating a request is byte-identical.

Fewer than two run ids in the path is a 400 `validation_error` with the
message `cohort requires at least two distinct runs`; a repeated run id is a
409 `conflict` with the message `cohort run ids must be distinct`; any run id
that does not exist is a 404 `not_found` with the message
`run <run_id> was not found`.

### The annotation table

Annotation matches the exact tuple `(CHROM, POS, REF, ALT)`. A variant with no
matching row is annotated `gene: null`, `consequence: "intergenic_variant"`,
`impact: "MODIFIER"`.

```
chr1  11856378  G  A    MTHFR  missense_variant       MODERATE
chr1  11856378  G  T    MTHFR  stop_gained            HIGH
chr7  55019017  G  GA   EGFR   frameshift_variant     HIGH
chr7  55019017  G  T    EGFR   stop_gained            HIGH
chr7  55273000  C  T    EGFR   synonymous_variant     LOW
chr11 5227002   T  C    HBB    missense_variant       MODERATE
chr12 25245350  C  A    KRAS   missense_variant       MODERATE
chr12 25245350  C  T    KRAS   missense_variant       MODERATE
chr12 25245351  C  T    KRAS   splice_donor_variant   HIGH
chr17 43093445  C  T    BRCA1  stop_gained            HIGH
chr17 43093446  A  G    BRCA1  missense_variant       MODERATE
chr17 7676154   C  T    TP53   missense_variant       MODERATE
chr19 11200138  C  T    LDLR   missense_variant       MODERATE
chrX  31137345  A  G    DMD    missense_variant       MODERATE
```

This is a small synthetic table with GRCh38-style chromosome names and real gene
symbols; the coordinates are illustrative and it is not a clinical annotation
source.

### Get run provenance

`GET /runs/run-1/provenance` returns `{"run_id","verified","chain_verified",
"reproduction_verified","head_sha256","steps":[...]}`, where a step looks like:

```json
{"step":1,"name":"ingest","params":{"parser_version":"variantrail-vcf-1"},
 "input_sha256":"45520c9b2b4798adafe3bdbaf137c567536527581c259dc2d92c391baf2a94ef",
 "output_sha256":"0f5d73dd61762cefc2df233d1e604de50d21df017a6f9e5608851a3c700c6103",
 "previous_sha256":null,
 "hash":"6de42193fc25e82afb9c32b0d6170629b69eb25268d978864f9162b27ed75530"}
```

Write `canonical(x)` for `json.dumps(x, ensure_ascii=False, separators=(",",":"),
sort_keys=True)` and `H(x)` for the lowercase SHA-256 hex digest of the UTF-8
bytes of `canonical(x)`. The four links are:

| step | name | params | input_sha256 | output_sha256 |
| --- | --- | --- | --- | --- |
| 1 | `ingest` | `{"parser_version":"variantrail-vcf-1"}` | SHA-256 of the VCF text | `H({"records":<parsed records>})` |
| 2 | `filter` | `min_dp`, `min_qual`, `pass_only`, plus the normalized `record_filter` when supplied | step 1 `output_sha256` | `H({"filter_counts":<record-reason keys>,"records":<kept records>})` |
| 3 | `annotate` | `annotation_table_sha256`, `genes`, `impacts` | step 2 `output_sha256` | `H({"filter_counts":<record keys plus allele_filtered>,"variants":<variants>})` |
| 4 | `summarize` | `{}` | step 3 `output_sha256` | `H(<statistics>)` |

Without a `record_filter` the step-2 filter counts hold the three original
record-reason keys and the step params hold only `min_dp`, `min_qual` and
`pass_only`, so every link hash is unchanged from the five-parameter contract.
With one, `rule_mismatch` is added (after `dp_below_min`) and the normalized
expression is part of the step-2 params; reproduction re-runs the saved
expression and the rule result, statistics, variants and all four step output
hashes feed the existing `verified` verdict.

Each link's own hash covers the previous link, so the trail is tamper-evident:

```
body = {step, name, params, input_sha256, output_sha256, previous_sha256}
hash = H(body)                 # previous_sha256 is null for step 1
```

`chain_verified` recomputes every link hash and checks that each
`previous_sha256` equals the previous link's `hash`. `reproduction_verified`
re-reads the stored sample VCF, re-executes all four steps with the stored
`params`, and compares the four `output_sha256` values, the statistics and the
variants. `verified` is the conjunction of both; the endpoint never rewrites the
stored trail.

### Step-level lineage manifest

The run entry point `variantrail.pipeline.execute(parsed, vcf_text, params,
lineage=False)` accepts an optional boolean `lineage` argument. It is off by
default, and an off run returns exactly the legacy result object. With
`lineage=True` the result gains one extra key, `lineage`, a self-contained
manifest that needs no logs or stored files to explain the result:

```json
{"input":{"sha256":"<sha256 of the VCF text>"},
 "params":{<complete parameter snapshot taken at call time>},
 "steps":[{"input_sha256":"<digest of the step input>",
           "name":"ingest",
           "output_sha256":"<digest of the step output>",
           "params":{<parameters the step actually used>},
           "previous_sha256":"0000...0000",
           "sha256":"<this record's digest>"},
          ...],
 "final_sha256":"<digest of the last record in the chain>"}
```

- `input.sha256` is the SHA-256 of the VCF text and equals the first record's
  `input_sha256`.
- `params` is a complete snapshot of the normalized parameters deep-copied via a
  JSON round trip before any step runs. Mutating the caller's parameter object
  after the call cannot change the returned manifest, and a parameter value that
  cannot be represented as JSON makes `execute` raise `TypeError` before any
  step executes.
- `steps` lists only the steps that actually executed, in real execution order.
  An empty input or a filter that removes every variant still records the steps
  that ran; a step skipped by a condition never appears as an executed record.
- The first record's `previous_sha256` is exactly 64 `0` characters; every
  later record carries the previous record's `sha256`, and `final_sha256`
  equals the last record's `sha256`.

Each record digest is the SHA-256 of the record body over every field except
`sha256` itself, rendered as canonical JSON: `json.dumps(body,
ensure_ascii=False, separators=(",", ":"), sort_keys=True)` encoded as UTF-8.
Object keys are sorted by Unicode code point, array order is preserved, and
numbers, booleans and `null` keep their JSON types. Because every input digest
is canonical, dictionaries whose keys are ordered differently but whose values
are equal hash identically.

The manifest is deterministic: identical input content, equivalent parameters
and identical step outputs produce a byte-identical manifest on repeated runs,
and changing any input, parameter or step output changes every record digest
from the first affected record through `final_sha256`.

### Verifying a lineage manifest

`variantrail.verify_lineage(manifest)` (also exported from the package root) is
the public verification entry point. It recomputes every record digest and the
chain linkage and returns `True` only when the whole manifest is consistent:
each `previous_sha256` chains correctly from the 64-zero head, each `sha256`
recomputes from the other fields, the first input digest matches `input.sha256`,
the top-level parameter snapshot matches the parameters carried by the steps,
and `final_sha256` closes the chain. It never mutates its argument.

Any tampered field, swapped or missing record, wrong chain-head value or
mismatching `final_sha256` returns `False`. A non-object argument, a missing
required field, or a digest that is not a 64-character lowercase hexadecimal
string raises `ValueError`.

## Errors

```json
{"error":{"code":"validation_error","message":"line 5: POS must be a positive integer"}}
```

Validation errors return 400, missing resources 404, conflicts 409. A missing
`Idempotency-Key` header is a `validation_error`, and reusing one key for a
different operation is a `conflict`.

## Tests

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```
