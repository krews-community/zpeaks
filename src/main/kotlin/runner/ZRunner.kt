package runner

import io.createSignalFile
import io.writeSkewSubPeaksBed
import io.writeStandardSubPeaksBed
import model.*
import step.*
import mu.KotlinLogging
import step.subpeaks.*
import java.nio.file.Path
import java.util.concurrent.ForkJoinPool

private val log = KotlinLogging.logger {}

data class ZRunConfig(
    val pileUpInputs: List<PileUpInput>,
    val chrFilter: Map<String, IntRange?>? = null,
    val signalOut: SignalOutput? = null,
    val peaksOut: Path? = null,
    val smoothing: Double = 50.0,
    val threshold: Double = 6.0,
    val fitMode: FitMode = FitMode.SKEW,
    val parallelism: Int? = null
)

abstract class ZRunner(private val name: String, protected val runConfig: ZRunConfig) {

    fun prepBams() = with(runConfig) { prepBams(pileUpInputs.map { it.bam }, chrFilter) }
    abstract fun pileUp(chr: String, chrLength: Int, range: IntRange): PileUp
    fun pdf(pileUp: PileUp) = with(runConfig) { pdf(pileUp, smoothing) }
    open fun peaks(pdf: PDF): List<Region> = with(runConfig) { callPeaks(pdf, threshold) }

    fun run() = with(runConfig) {
        if (parallelism != null) {
            System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", parallelism.toString())
        }
        log.info { "ZPeaks run started with parallelism = ${ForkJoinPool.commonPool().parallelism}" }

        log.info { "Running $name ZPeaks operation on ${pileUpInputs.size} alignment files..." }
        val chrsWithBounds = prepBams()
        for ((chr, bounds) in chrsWithBounds) {
            val pileUp = pileUp(chr, bounds.length, bounds.range)
            if (signalOut != null && signalOut.type == SignalOutputType.RAW) {
                createSignalFile(signalOut.path, signalOut.format, chr, pileUp)
                if (peaksOut == null) continue
            }

            val pdf = pdf(pileUp)
            if (0.0 == pdf.background.average || 0.0 == pdf.background.stdDev) {
                log.warn { "One or more background parameters for chromosome $chr was zero. skipping..." }
                continue
            }

            if (signalOut != null && signalOut.type == SignalOutputType.SMOOTHED) {
                createSignalFile(signalOut.path, signalOut.format, chr, pdf, signalOut.signalResolution)
            }

            if (peaksOut == null) continue
            val peaks = peaks(pdf)
            if (fitMode == FitMode.SKEW) {
                val skewSubPeaks = SkewFitter.fit(chr, peaks, pdf)
                writeSkewSubPeaksBed(peaksOut, chr, skewSubPeaks)
            } else {
                val subPeaks = StandardFitter.fit(chr, peaks, pdf)
                writeStandardSubPeaksBed(peaksOut, chr, subPeaks)
            }
        }
        log.info { "$name ZPeaks run complete!" }
    }
}