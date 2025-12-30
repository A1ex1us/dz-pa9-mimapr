import java.io.File
import kotlin.math.abs
import kotlin.math.exp
import kotlin.math.sqrt

fun gauss(a: Array<DoubleArray>, b: DoubleArray): DoubleArray {
    if (a.size != b.size) throw MatrixDimensionException()

    val x = DoubleArray(b.size) { 0.0 }

    for (j in 0..<a.size - 1) {
        val ajj = a[j][j]
        for (i in j + 1..<a.size) {
            val coeff = a[i][j] / ajj
            for (k in j..<a.size) {
                a[i][k] -= coeff * a[j][k]
            }
            b[i] -= coeff * b[j]
        }
    }
    for (i in a.size - 1 downTo 0) {
        var sum = 0.0
        for (ind in a.size - 1 downTo i + 1) {
            sum += a[i][ind] * x[ind]
        }
        x[i] = (b[i] - sum) / a[i][i]
    }
    return x
}

fun createJacobiMatrix(deltaT: Double, uCb1: Double, uCb2: Double): Array<DoubleArray> {
    val jacobi = Array<DoubleArray>(Data.N) { DoubleArray(Data.N) { 0.0 } }

    with(Components) {
        val a = -1.0 / R_U - I_T / MFT * exp(uCb1 / MFT)
        val b = 1.0 / R_U + I_T / MFT * exp(uCb2 / MFT)

        jacobi[0][0] = 1.0
        jacobi[0][5] = -1.0 / deltaT
        jacobi[1][1] = 1.0
        jacobi[1][6] = -1.0 / deltaT
        jacobi[2][2] = 1.0
        jacobi[2][7] = -1.0 / deltaT
        jacobi[3][3] = 1.0
        jacobi[3][8] = -1.0 / deltaT
        jacobi[4][4] = 1.0
        jacobi[4][9] = -1.0 / deltaT

        jacobi[5][5] = 1.0
        jacobi[5][10] = -1.0
        jacobi[5][11] = 1.0

        jacobi[6][6] = 1.0
        jacobi[6][13] = -1.0

        jacobi[7][7] = 1.0
        jacobi[7][11] = 1.0
        jacobi[7][15] = -1.0

        jacobi[8][8] = 1.0
        jacobi[8][12] = -1.0
        jacobi[8][13] = 1.0

        jacobi[9][4] = L
        jacobi[9][13] = -1.0
        jacobi[9][14] = 1.0

        jacobi[10][0] = C1
        jacobi[10][16] = 1.0

        jacobi[11][0] = -C1
        jacobi[11][2] = -C_B
        jacobi[11][7] = a
        jacobi[11][11] = 1.0 / R_B
        jacobi[11][12] = -1.0 / R_B

        jacobi[12][3] = C_B
        jacobi[12][8] = b
        jacobi[12][11] = -1.0 / R_B
        jacobi[12][12] = 1.0 / R_B

        jacobi[13][1] = C2
        jacobi[13][3] = -C_B
        jacobi[13][8] = -b
        jacobi[13][9] = 1.0
        jacobi[13][13] = 1.0 / R

        jacobi[14][9] = -1.0
        jacobi[14][17] = 1.0

        jacobi[15][2] = C_B
        jacobi[15][7] = -a
        jacobi[15][15] = 1.0 / R_B

        jacobi[16][10] = -1.0
        jacobi[17][14] = -1.0
    }
    return jacobi
}

fun createVector(
    td: TimeDemon,
    pvApprox: PhaseVariables,
    prevStateVariables: PrevStateVariables
): DoubleArray {
    // распаковка объекта приближений фазовых переменных в отдельные значения
    val (dUC1dt, dUC2dt, dUCB1dt, dUCB2dt, dILdt, uC1, uC2,
        uCB1, uCB2, iL, phi1, phi2, phi3, phi4, phi5, phi6, iE1, iE2) = pvApprox

    val (uC1Prev, uC2Prev, uCB1Prev, uCB2Prev, iLPrev) = prevStateVariables

    val vector = DoubleArray(Data.N) { 0.0 }

    val deltaT = td.deltaT

    with(Components) {
        vector[0] = dUC1dt - (uC1 - uC1Prev) / deltaT
        vector[1] = dUC2dt - (uC2 - uC2Prev) / deltaT
        vector[2] = dUCB1dt - (uCB1 - uCB1Prev) / deltaT
        vector[3] = dUCB2dt - (uCB2 - uCB2Prev) / deltaT
        vector[4] = dILdt - (iL - iLPrev) / deltaT

        vector[5] = uC1 - (phi1 - phi2)
        vector[6] = uC2 - phi4
        vector[7] = uCB1 - (phi6 - phi2)
        vector[8] = uCB2 - (phi3 - phi4)
        vector[9] = L * dILdt - (phi4 - phi5)

        vector[10] = iE1 + C1 * dUC1dt
        vector[11] = -C1 * dUC1dt - C_B * dUCB1dt - uCB1 / R_U - I_T * (exp(uCB1 / MFT) - 1.0) + (phi2 - phi3) / R_B
        vector[12] = -(phi2 - phi3) / R_B + C_B * dUCB2dt + uCB2 / R_U + I_T * (exp(uCB2 / MFT) - 1.0)
        vector[13] = -(C_B * dUCB2dt + uCB2 / R_U + I_T * (exp(uCB2 / MFT) - 1.0)) + C2 * dUC2dt + phi4 / R + iL
        vector[14] = -iL + iE2
        vector[15] = phi6 / R_B + C_B * dUCB1dt + uCB1 / R_U + I_T * (exp(uCB1 / MFT) - 1.0)

        vector[16] = getCurrentE1(td.currT) - phi1
        vector[17] = E2 - phi5
    }

    return vector
}

fun newtonMethod(
    td: TimeDemon, initApprox: PhaseVariables,
    prevStateVariables: PrevStateVariables
): NewtonMethodResults {
    var n = 0

    var currApprox = initApprox.copy()
    while (n < Data.NEWTON_STEPS) {

        val jacobi = createJacobiMatrix(
            deltaT = td.deltaT,
            uCb1 = initApprox.uCb1,
            uCb2 = initApprox.uCb2
        )
        val vectorForNewton = createVector(
            td = td,
            pvApprox = currApprox,
            prevStateVariables = prevStateVariables
        ).apply {
            for (i in this.indices) {
                this[i] *= -1.0
            }
        }

        val gaussResults = gauss(jacobi, vectorForNewton)
        val deltas = PhaseVariables.convertFromArrayToPV(gaussResults)
        var newApprox = currApprox + deltas
        currApprox = newApprox

        if (calculateVectorNorm(deltas) < Data.EPSILON) break

        n++
    }
    if (n >= 7) {
        return NewtonMethodResults(false)
    }

    return NewtonMethodResults(true, currApprox)
}

private fun findMaxValue(values: DoubleArray): Double {
    var max = 0.0
    values.forEach { component ->
        if (component > max) max = component
    }
    return max
}

private fun calculateVectorNorm(pv: PhaseVariables): Double {
    val array = pv.convertFromPVToArray()
    val result = array.fold(0.0) { acc: Double, value -> acc + value * value }
    return sqrt(result)
}

fun predictPhaseVariables(
    pvPrev: PhaseVariables,
    pvPrevPrev: PhaseVariables,
): PhaseVariables {
    val prevPrevArray = pvPrevPrev.convertFromPVToArray()
    val newValues = pvPrev
        .convertFromPVToArray()
        .mapIndexed { i, value -> 2.0 * value - prevPrevArray[i] }
        .toDoubleArray()
    return PhaseVariables.convertFromArrayToPV(newValues)
}

fun calculateDeltaT(
    td: TimeDemon, prevDeltaT: Double,
    pv: PhaseVariables, pvPrev: PhaseVariables, pvPrevPrev: PhaseVariables
): TimeDemon {
    val currT = td.currT
    val deltaT = td.deltaT

    val secondRemainder: (Double, Double, Double) -> Double = { curr, prev, prevPrev ->
        abs((curr * prevDeltaT - prev * (deltaT + prevDeltaT) + prevPrev * deltaT) / (2 * prevDeltaT))
    }

    val delta: Double = findMaxValue(
        PhaseVariables(
            phi1 = secondRemainder(pv.phi1, pvPrev.phi1, pvPrevPrev.phi1),
            phi2 = secondRemainder(pv.phi2, pvPrev.phi2, pvPrevPrev.phi2),
            phi3 = secondRemainder(pv.phi3, pvPrev.phi3, pvPrevPrev.phi3),
            phi4 = secondRemainder(pv.phi4, pvPrev.phi4, pvPrevPrev.phi4),
            phi5 = secondRemainder(pv.phi5, pvPrev.phi5, pvPrevPrev.phi5),
            phi6 = secondRemainder(pv.phi6, pvPrev.phi6, pvPrevPrev.phi6),
        ).convertFromPVToArray()
    )

    return if (delta > Data.DELTA_1) {
        TimeDemon(currT = currT, deltaT = deltaT / 2, success = false)
    } else {
        if (delta > Data.DELTA_2) {
            TimeDemon(currT = currT + deltaT, deltaT = deltaT, success = true)
        } else {
            val newDeltaT = if (deltaT > Data.MAXIMAL_STEP) Data.MAXIMAL_STEP else deltaT * 2
            TimeDemon(currT = currT + deltaT, deltaT = newDeltaT, success = true)
        }
    }
}


fun deltaTReduction(timeDemon: TimeDemon) =
    timeDemon.copy(deltaT = timeDemon.deltaT / 2.0)

fun printToFile(name: String, values: DoubleArray) {
    val file = File(name)
    for (i in values.indices) {
        file.appendText("${values[i]}\n")
    }
}

fun clearFiles() {
    val files = listOf(
        File(FileData.PHI1_FILE), File(FileData.PHI2_FILE),
        File(FileData.PHI4_FILE), File(FileData.PHI5_FILE), File(FileData.T_FILE),
        // Новые файлы
        File(FileData.UC1_FILE),
        File(FileData.UC2_FILE),
        File(FileData.IL_FILE)
    )
    files.forEach { file -> if (file.exists()) file.delete() }
}
