import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.*
import io.*
import model.*
import mu.KotlinLogging
import step.*
import step.subpeaks.*
import util.*
import java.nio.file.Path
import java.util.*

private val log = KotlinLogging.logger {}

class ZPeaks : CliktCommand() {
    val samIn: Path by option("-samIn", help="Input Sam or Bam alignment file").path().required()
    val pileUpOut: Path? by option("-pileUp", help="Output Pile-Up file").path()
    val pileUpOutFormat: PileUpOutputFormat? by option("")
        .choice(PileUpOutputFormat.values().associateBy { it.lowerHyphenName })
        .required()
    val peaksOut: Path? by option("-peaksOut", help="Output Peaks bed file").path()
    val subPeaksOut: Path? by option("-subPeaksOut", help="Output Sub-Peaks bed file").path()
    val strand: Strand by option("-strand", help="Strand to count during pile-up")
        .choice(Strand.values().associateBy { it.lowerHyphenName })
        .default(Strand.BOTH)
    val atacMode: Boolean by option("-atacMode",
        help="Whether or not this alignment is for an ATAC experiment. Algorithm will use special offsets if true.")
        .flag()
    val pileUpAlgorithm: PileUpAlgorithm by option("-pileUpAlgorithm",
        help="Algorithm used to select values during pile up.")
        .choice(PileUpAlgorithm.values().associateBy { it.lowerHyphenName })
        .default(PileUpAlgorithm.START)
    val smoothing: Double by option("-smoothingFactor",
        help="Smoothing Factor for calculating PDF for Pile Up data during Peaks step.").double()
        .default(50.0)
    val noAmplitude: Boolean by option("-noAmplitude",
        help="Used during PDF calculation during Peaks step").flag() //TODO add more to this help description
    val threshold: Double by option("-threshold", help="Threshold used during peak calling").double()
        .default(6.0)

    override fun run() {
        val pileUpOutput = if (pileUpOut != null) PileUpOutput(pileUpOut!!, pileUpOutFormat!!) else null
        val pileUps = runPileUp(samIn, strand, atacMode, pileUpAlgorithm, pileUpOutput)
        runPeaks(pileUps, smoothing, noAmplitude, threshold, peaksOut, subPeaksOut)
    }
}

fun main(args: Array<String>) = ZPeaks().main(args)

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
        subPeaks += fitSkew(peakValues, region.start).flatMap { fit ->
            fit.subPeaks
        }
        val elapsed = (System.currentTimeMillis() - start) / 1000.0
        log.info { "Job took $elapsed seconds" }
    }
    subPeaks.sortBy { it.region.start }
    return subPeaks
}