import com.github.ajalt.clikt.core.CliktCommand
import io.*
import model.*
import mu.KotlinLogging
import step.*
import step.subpeaks.*
import util.runParallel
import java.nio.file.Path
import java.util.*

private val log = KotlinLogging.logger {}

class ZPeaks : CliktCommand() {


    override fun run() {

    }
}

fun main(args: Array<String>) = ZPeaks().main(args)

fun runPeaks(pileUps: Map<String, PileUp>, smoothing: Double, noAmplitude: Boolean, threshold: Double,
             peaksOut: Path? = null, subPeaksOut: Path? = null, onRange: IntRange? = null) {
    for ((chr, pileUp) in pileUps) {
        log.info { "Calculating PDF for chromosome $chr pileup data..." }
        val pdf = pdf(chr, pileUp, smoothing, noAmplitude, onRange)
        log.info { "Chromosome $chr PDF completed with background ${pdf.background}" }

        if (0.0 == pdf.background.average || 0.0 == pdf.background.stdDev) {
            log.warn { "One or more background parameters for chromosome $chr was zero. skipping..." }
            continue
        }

        val peaks = callPeaks(pdf, threshold)
        log.info { "Peaks called for threshold $threshold. Peaks found: ${peaks.size}" }

        if (peaksOut != null) writePeaksBed(peaksOut, chr, peaks)
        if (subPeaksOut == null) continue

        val subPeaks = runSkewSubPeaks(pdf, peaks as MutableList, 6)
        writeSkewSubPeaksBed(subPeaksOut, chr, subPeaks)
    }
}

fun runSkewSubPeaks(pdf: PDF, peaks: List<Peak>, parallelism: Int): List<SkewSubPeak> {
    (peaks as MutableList).sortByDescending { it.region.end - it.region.start }
    val subPeaks = Collections.synchronizedList(mutableListOf<SkewSubPeak>())
    runParallel("Skew Sub-Peaks", peaks, parallelism) { peak ->
        val start = System.currentTimeMillis()
        val region = peak.region
        val peakValues = (region.start..region.end).map { pdf[it] }
        subPeaks += fitSkew(peakValues, region.start)
        val elapsed = (System.currentTimeMillis() - start) / 1000.0
        log.info { "Job took $elapsed seconds" }
    }
    subPeaks.sortBy { it.region.start }
    return subPeaks
}

data class PileUpOutput(val path: Path, val format: PileUpOutputFormat)
enum class PileUpOutputFormat { WIG, BIG_WIG, BED_GRAPH }

fun runPileUp(samIn: Path, strand: Strand, atacMode: Boolean, pileUpAlgorithm: PileUpAlgorithm,
              pileUpOut: PileUpOutput? = null): Map<String, PileUp> {
    log.info { "Performing pile-up for sam file $samIn with strand=$strand, atacMode=$atacMode, pileUpAlgorithm=$pileUpAlgorithm" }
    val pileUps = pileUpSam(samIn, strand, atacMode, pileUpAlgorithm)

    val pileUpSummary = pileUps.entries.joinToString("\n") {
        "chromosome ${it.key}: len=${it.value.chromosomeLength} sum=${it.value.sum}"
    }
    log.info { "Pile-up complete with results: \n$pileUpSummary" }

    if (pileUpOut != null) when(pileUpOut.format) {
        PileUpOutputFormat.WIG -> { writeWig(pileUpOut.path, pileUps) }
        PileUpOutputFormat.BIG_WIG -> { writeBigWig(pileUpOut.path, pileUps) }
        PileUpOutputFormat.BED_GRAPH -> { writeBedGraph(pileUpOut.path, pileUps) }
    }

    return pileUps
}
