import io.createSignalFile
import io.writeSkewSubPeaksBed
import io.writeStandardSubPeaksBed
import model.*
import mu.KotlinLogging
import step.*
import step.subpeaks.SkewFitter
import step.subpeaks.StandardFitter
import java.nio.file.Path

private val log = KotlinLogging.logger {}

fun runBottomUp(pileUpInputs: List<PileUpInput>, chrFilter: List<String>?, signalOut: SignalOutput? = null,
                smoothing: Double, normalizePDF: Boolean, threshold: Double, peaksOut: Path?,
                fitMode: FitMode) {
    log.info { "Running bottom-up peaks operation on ${pileUpInputs.size} alignment files..." }
    val chrLengths = prepBams(pileUpInputs.map { it.bam }, chrFilter)
    for ((chr, length) in chrLengths) {
        log.info { "Creating aggregate pile-up for $chr..." }
        val aggregatePileUpData = DoubleArray(length)
        var sum = 0.0
        for (pileUpInput in pileUpInputs) {
            val pileUp = runPileUp(pileUpInput.bam, chr, length, pileUpInput.options)
            for (index in 0 until length) {
                val value = pileUp.scaledValue(index)
                aggregatePileUpData[index] += value
                sum += value
            }
        }

        log.info { "Aggregate pile-up for $chr complete!" }
        val aggregatePileUp = PileUp(aggregatePileUpData, chr, length, sum)

        if (signalOut != null && signalOut.type == SignalOutputType.RAW) {
            createSignalFile(signalOut.path, signalOut.format, chr, aggregatePileUp)
        }

        log.info { "Smoothing aggregate pile-up for $chr..." }
        val pdf = pdf(aggregatePileUp, smoothing, normalizePDF)
        log.info { "Aggregate $chr PDF completed with background ${pdf.background}" }

        if (0.0 == pdf.background.average || 0.0 == pdf.background.stdDev) {
            log.warn { "One or more background parameters for chromosome $chr was zero. skipping..." }
            continue
        }

        if (signalOut != null && signalOut.type == SignalOutputType.SMOOTHED) {
            createSignalFile(signalOut.path, signalOut.format, chr, pdf, signalOut.signalResolution)
        }

        log.info { "Calling peaks on smoothed aggregate $chr data..." }
        val peaks = callPeaks(pdf, threshold)
        log.info { "Peaks call complete! ${peaks.size} peaks found." }

        if (peaksOut == null) return
        if (fitMode == FitMode.SKEW) {
            val skewSubPeaks = SkewFitter.fit(chr, peaks, pdf)
            writeSkewSubPeaksBed(peaksOut, chr, skewSubPeaks)
        } else {
            val subPeaks = StandardFitter.fit(chr, peaks, pdf)
            writeStandardSubPeaksBed(peaksOut, chr, subPeaks)
        }
    }
    log.info { "Bottom-up peaks run complete!" }
}