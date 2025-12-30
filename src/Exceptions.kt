interface MimaprException {
    val error: String
}

data class MatrixDimensionException(
    override val error: String = "Invalid size of matrix"
) : RuntimeException(), MimaprException

data class TimeStepException(
    override val error: String = "Step is now less than the minimum possible"
) : RuntimeException(), MimaprException

data class ConversionException(
    override val error: String = "Gauss results don't match the required size"
): RuntimeException(), MimaprException