import io.createSignalFile
import model.*
import mu.KotlinLogging
import step.*
import java.nio.file.Path

private val log = KotlinLogging.logger {}

fun runTopDown(pileUpInputs: List<PileUpInput>, chrFilter: List<String>?, signalOut: SignalOutput? = null,
                smoothing: Double, normalizePDF: Boolean, threshold: Double, peaksOut: Path?,
                fitMode: FitMode) {
    log.info { "Running top-down peaks operation on ${pileUpInputs.size} alignment files..." }
    val chrLengths = prepBams(pileUpInputs.map { it.bam }, chrFilter)
    for ((chr, length) in chrLengths) {
        log.info { "Creating peaks-only aggregate track..." }
        // Array containing all aggregated data
        val aggregateChrData = DoubleArray(length)
        // Array containing the number of sources the data came from at each BP
        val aggregateSources = IntArray(length)
        val allPeaks = mutableListOf<List<Region>>()

        for (pileUpInput in pileUpInputs) {
            val pileUp = runPileUp(pileUpInput.bam, chr, length, pileUpInput.options)
            val pdf = pdf(pileUp, smoothing, normalizePDF)
            val peaks = callPeaks(pdf, threshold)
            allPeaks += peaks
            // Add raw pile-up data to aggregate from peak regions only
            for (peak in peaks) {
                for (bp in peak.start .. peak.end) {
                    aggregateChrData[bp] += pileUp.scaledValue(bp)
                    aggregateSources[bp]++
                }
            }
        }

        // Scale the aggregate peaks data based on the number of sources each bp came from
        var aggregateDataSum = 0.0
        for (bp in 0 until aggregateChrData.size) {
            aggregateChrData[bp] /= aggregateSources[bp].toDouble()
            aggregateDataSum += aggregateChrData[bp]
        }

        log.info { "Peaks-only aggregate Created! Smoothing aggregate pile-up for $chr..." }
        val aggregatePileUp = PileUp(aggregateChrData, chr, length, aggregateDataSum)
        val pdf = pdf(aggregatePileUp, smoothing, normalizePDF)
        log.info { "Aggregate $chr PDF completed with background ${pdf.background}" }

        if (0.0 == pdf.background.average || 0.0 == pdf.background.stdDev) {
            log.warn { "One or more background parameters for chromosome $chr was zero. skipping..." }
            continue
        }

        if (signalOut != null && signalOut.type == SignalOutputType.SMOOTHED) {
            createSignalFile(signalOut.path, signalOut.format, chr, pdf, signalOut.signalResolution)
        }

        val mergedPeaks = mergePeaks(allPeaks)

    }
    log.info { "Top-down peaks run complete!" }
}