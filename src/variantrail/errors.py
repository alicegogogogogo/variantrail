class VariantRailError(Exception):
    code = "internal_error"
    status = 500


class ValidationError(VariantRailError):
    code = "validation_error"
    status = 400


class NotFoundError(VariantRailError):
    code = "not_found"
    status = 404


class ConflictError(VariantRailError):
    code = "conflict"
    status = 409
