fun main() {
    clearFiles()

    var timeDemon = TimeDemon(
        currT = Data.TIME_START,
        deltaT = Data.START_DELTA_TIME
    )
    var prevDeltaT = timeDemon.deltaT

    var initialApproximation = PhaseVariables()
    var pvPrev = PhaseVariables()
    var pvPrevPrev = PhaseVariables()

    val results = ResultLists()

    val addResults: (PhaseVariables, Double) -> Unit = { pv, time ->
        results.phi1List.add(pv.phi1)
        results.phi2List.add(pv.phi2)
        results.phi4List.add(pv.phi4)
        results.phi5List.add(pv.phi5)
        results.timeList.add(time)
        // Добавление новых данных
        results.uC1List.add(pv.uC1)   // Напряжение на C1
        results.uC2List.add(pv.uC2)   // Напряжение на C2
        results.iLList.add(pv.iL)     // Ток через индуктивность
    }
    addResults(initialApproximation, timeDemon.currT)

    var prevStateVariables = PrevStateVariables(
        uC1 = 0.0,
        uC2 = 0.0,
        uCb1 = 0.0,
        uCb2 = 0.0,
        iL = 0.0
    )

    var iteration = 0
    while (timeDemon.currT < Data.TIME_END) {
        var newtonMethodResults = newtonMethod(
            td = timeDemon,
            initApprox = initialApproximation,
            prevStateVariables = prevStateVariables
        )
        iteration++

        when (newtonMethodResults.isSuccessful) {
            true -> {
                val prevTimeDemon = timeDemon.copy()
                timeDemon = calculateDeltaT(
                    td = timeDemon,
                    prevDeltaT = prevDeltaT,
                    pv = newtonMethodResults.phaseVariables,
                    pvPrev = pvPrev,
                    pvPrevPrev = pvPrevPrev
                )

                if (!timeDemon.success) {
                    continue
                }

                initialApproximation = predictPhaseVariables(newtonMethodResults.phaseVariables, pvPrev)

                with(newtonMethodResults.phaseVariables) {
                    prevStateVariables = PrevStateVariables(
                        uC1 = uC1,
                        uC2 = uC2,
                        uCb1 = uCb1,
                        uCb2 = uCb2,
                        iL = iL
                    )
                }

                pvPrevPrev = pvPrev.copy()
                pvPrev = newtonMethodResults.phaseVariables.copy()
                if (timeDemon.currT >= Data.TIME_END - timeDemon.deltaT || iteration % 1000 == 0) {
                    println("$iteration, t=${timeDemon.currT}")
                    addResults(newtonMethodResults.phaseVariables, timeDemon.currT)
                }
            }

            false -> {
                timeDemon = deltaTReduction(timeDemon)

                if (timeDemon.deltaT / 2.0 < Data.MINIMAL_STEP) throw TimeStepException()
            }
        }
    }

    with(results) {
        printToFile(FileData.PHI1_FILE, phi1List.toDoubleArray())
        printToFile(FileData.PHI2_FILE, phi2List.toDoubleArray())
        printToFile(FileData.PHI4_FILE, phi4List.toDoubleArray())
        printToFile(FileData.PHI5_FILE, phi5List.toDoubleArray())
        printToFile(FileData.T_FILE, timeList.toDoubleArray())
        // Сохранение новых данных
        printToFile(FileData.UC1_FILE, uC1List.toDoubleArray())   // C1
        printToFile(FileData.UC2_FILE, uC2List.toDoubleArray())   // C2
        printToFile(FileData.IL_FILE, iLList.toDoubleArray())     // Индуктивность
    }
}
