import kotlin.math.PI
import kotlin.math.sin

data object Data {
    const val START_DELTA_TIME = 1e-11
    const val TIME_START = 0.0
    const val TIME_END = 1e-3
    const val MINIMAL_STEP = 1e-18
    const val MAXIMAL_STEP = 1e-5

    const val N = 18
    const val NEWTON_STEPS = 7

    const val EPSILON = 5e-2
    const val DELTA_1 = 1e-2
    const val DELTA_2 = 1e-3
}

data object Components {
    const val L = 2.53e-4
    const val C1 = 1e-6
    const val C2 = 1e-6
    const val I_T = 1e-12
    const val C_B = 2e-12
    const val MFT = 0.026
    const val R_U = 1_000_000.0
    const val R_B = 20.0
    const val R = 1_000.0
    val getCurrentE1: (Double) -> Double = { t -> 10.0 * sin(2.0 * PI / 1e-4 * t) }
    const val E2 = 5.0
}

data object FileData {
    const val PHI1_FILE = "phi1.txt"
    const val PHI2_FILE = "phi2.txt"
    const val PHI4_FILE = "phi4.txt"
    const val PHI5_FILE = "phi5.txt"
    const val T_FILE = "t.txt"
    // Новые файлы для сохранения данных элементов
    const val UC1_FILE = "uC1.txt"   // Напряжение на конденсаторе C1
    const val UC2_FILE = "uC2.txt"   // Напряжение на конденсаторе C2
    const val IL_FILE = "iL.txt"     // Ток через индуктивность L
}

data class PhaseVariables(
    val dUC1dt: Double = 0.0,
    val dUC2dt: Double = 0.0,
    val dUCB1dt: Double = 0.0,
    val dUCB2dt: Double = 0.0,
    val dILdt: Double = 0.0,

    val uC1: Double = 0.0,
    val uC2: Double = 0.0,
    val uCb1: Double = 0.0,
    val uCb2: Double = 0.0,
    val iL: Double = 0.0,

    val phi1: Double = 0.0,
    val phi2: Double = 0.0,
    val phi3: Double = 0.0,
    val phi4: Double = 0.0,
    val phi5: Double = 5.0,
    val phi6: Double = 0.0,

    val iE1: Double = 0.0,
    val iE2: Double = 0.0,
) {
    operator fun plus(other: PhaseVariables): PhaseVariables = PhaseVariables(
        dUC1dt = this.dUC1dt + other.dUC1dt,
        dUC2dt = this.dUC2dt + other.dUC2dt,
        dUCB1dt = this.dUCB1dt + other.dUCB1dt,
        dUCB2dt = this.dUCB2dt + other.dUCB2dt,
        dILdt = this.dILdt + other.dILdt,

        uC1 = this.uC1 + other.uC1,
        uC2 = this.uC2 + other.uC2,
        uCb1 = this.uCb1 + other.uCb1,
        uCb2 = this.uCb2 + other.uCb2,
        iL = this.iL + other.iL,

        phi1 = this.phi1 + other.phi1,
        phi2 = this.phi2 + other.phi2,
        phi3 = this.phi3 + other.phi3,
        phi4 = this.phi4 + other.phi4,
        phi5 = this.phi5 + other.phi5,
        phi6 = this.phi6 + other.phi6,

        iE1 = this.iE1 + other.iE1,
        iE2 = this.iE2 + other.iE2,
    )

    fun convertFromPVToArray(): DoubleArray = doubleArrayOf(
        this.dUC1dt,
        this.dUC2dt,
        this.dUCB1dt,
        this.dUCB2dt,
        this.dILdt,

        this.uC1,
        this.uC2,
        this.uCb1,
        this.uCb2,
        this.iL,

        this.phi1,
        this.phi2,
        this.phi3,
        this.phi4,
        this.phi5,
        this.phi6,

        this.iE1,
        this.iE2
    )

    companion object {
        fun convertFromArrayToPV(results: DoubleArray): PhaseVariables {
            if (results.size != Data.N) throw ConversionException()
            return PhaseVariables(
                dUC1dt = results[0],
                dUC2dt = results[1],
                dUCB1dt = results[2],
                dUCB2dt = results[3],
                dILdt = results[4],
                uC1 = results[5],
                uC2 = results[6],
                uCb1 = results[7],
                uCb2 = results[8],
                iL = results[9],
                phi1 = results[10],
                phi2 = results[11],
                phi3 = results[12],
                phi4 = results[13],
                phi5 = results[14],
                phi6 = results[15],
                iE1 = results[16],
                iE2 = results[17]
            )
        }
    }
}

data class PrevStateVariables(
    val uC1: Double,
    val uC2: Double,
    val uCb1: Double,
    val uCb2: Double,
    val iL: Double
)

data class NewtonMethodResults(
    val isSuccessful: Boolean,
    val phaseVariables: PhaseVariables = PhaseVariables()
)

data class TimeDemon(
    val currT: Double,
    val deltaT: Double,
    val success: Boolean = true
)

data class ResultLists(
    val phi1List: MutableList<Double> = mutableListOf(),
    val phi2List: MutableList<Double> = mutableListOf(),
    val phi4List: MutableList<Double> = mutableListOf(),
    val phi5List: MutableList<Double> = mutableListOf(),
    val timeList: MutableList<Double> = mutableListOf(),
    // Новые поля для элементов
    val uC1List: MutableList<Double> = mutableListOf(),
    val uC2List: MutableList<Double> = mutableListOf(),
    val iLList: MutableList<Double> = mutableListOf()
)