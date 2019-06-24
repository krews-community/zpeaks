import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.*
import io.*
import model.*
import mu.KotlinLogging
import step.*
import util.*
import java.nio.file.Path
import java.util.concurrent.ForkJoinPool


private val log = KotlinLogging.logger {}

class ZPeaksCommand(val run: (ZPeaksRunConfig) -> Unit = ::run): CliktCommand() {
    private val bamsIn: List<Path> by option("-bamIn", help="Input Bam alignment file")
        .path(exists = true)
        .multiple()
        .validate { require(it.isNotEmpty()) { "At least one path must be given" } }
    private val signalOutPath: Path? by option("-signalOut", help="Output Pile-Up file").path()
    private val signalOutType: SignalOutputType by option("-signalOutType")
        .choice(SignalOutputType.values().associateBy { it.lowerHyphenName })
        .default(SignalOutputType.SMOOTHED)
    private val signalOutFormat: SignalOutputFormat by option("-signalOutFormat")
        .choice(SignalOutputFormat.values().associateBy { it.lowerHyphenName })
        .default(SignalOutputFormat.BED_GRAPH)
    private val peaksOut: Path? by option("-peaksOut", help="Output Peaks bed file").path()
    private val strands: List<Strand> by option("-strand", help="Strand to count during pile-up")
        .choice(Strand.values().associateBy { it.lowerHyphenName })
        .multiple(listOf(Strand.BOTH))
    private val forwardShifts by option("-forwardShift", help="During pile-up, shift the " +
            "forward strand by this amount. Can be positive or negative. Default 0.")
        .int()
        .multiple(listOf(0))
    private val reverseShifts by option("-reverseShift", help="During pile-up, shift the " +
            "reverse strand by this amount. Can be positive or negative. Default 0.")
        .int()
        .multiple(listOf(0))
    private val pileUpAlgorithms: List<PileUpAlgorithm> by option("-pileUpAlgorithm",
        help="Algorithm used to select values during pile up.")
        .choice(PileUpAlgorithm.values().associateBy { it.lowerHyphenName })
        .multiple(listOf(PileUpAlgorithm.START))
    private val signalResolution: Int by option("-signalResolution",
        help="Number of decimal places to keep in outputted signal values")
        .int().default(1)
    private val smoothing: Double by option("-smoothingFactor",
        help="Smoothing Factor for calculating PDF for Pile Up data during Peaks step.").double()
        .default(50.0)
    private val normalizePDF: Boolean by option("-normalizePdf",
        help="If true, PDF values are scaled to the natural gaussian amplitude").flag()
    private val threshold: Double by option("-threshold", help="Threshold used during peak calling").double()
        .default(6.0)
    private val fitMode: FitMode by option("-fitMode", help="Sub-Peak Fitting Modes.")
        .choice(FitMode.values().associateBy { it.lowerHyphenName })
        .default(FitMode.SKEW)
    private val parallelism: Int? by option("-parallelism", help="Number of threads to use for parallel parts. " +
            "Defaults to number of cores on machine. Parallelism is NOT per-chromosome. Any amount is valid.").int()

    override fun run() {
        if (signalOutPath == null && peaksOut == null)
            throw Exception("One of the following must be set: -signalOut or -peaksOut")
        val signalOut =
            if (signalOutPath != null) SignalOutput(signalOutPath!!, signalOutType, signalOutFormat, signalResolution)
            else null

        val pileUpInputs = mutableListOf<PileUpInput>()
        for ((samIndex, sam) in bamsIn.withIndex()) {
            val strand = strands.getOrElse(samIndex) { strands.last() }
            val pileUpAlgorithm = pileUpAlgorithms.getOrElse(samIndex) { pileUpAlgorithms.last() }
            val forwardShift = forwardShifts.getOrElse(samIndex) { forwardShifts.last() }
            val reverseShift = reverseShifts.getOrElse(samIndex) { reverseShifts.last() }
            val pileUpOptions = PileUpOptions(strand, pileUpAlgorithm, forwardShift, reverseShift)
            pileUpInputs += PileUpInput(sam, pileUpOptions)
        }

        run(ZPeaksRunConfig(pileUpInputs, signalOut, peaksOut, smoothing, normalizePDF, threshold, fitMode, parallelism))
    }
}

enum class FitMode { SKEW, STANDARD }
data class SignalOutput(
    val path: Path,
    val type: SignalOutputType,
    val format: SignalOutputFormat,
    val signalResolution: Int = 1
)
enum class SignalOutputType { RAW, SMOOTHED }

data class ZPeaksRunConfig(
    val pileUpInputs: List<PileUpInput>,
    val signalOut: SignalOutput?,
    val peaksOut: Path?,
    val smoothing: Double,
    val normalizePDF: Boolean,
    val threshold: Double,
    val fitMode: FitMode = FitMode.SKEW,
    val parallelism: Int? = null
)

fun run(config: ZPeaksRunConfig) = with(config) {
    if (parallelism != null) {
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", parallelism.toString())
    }
    log.info { "ZPeaks run started with parallelism = ${ForkJoinPool.commonPool().parallelism}" }

    /*
    // Run peaks on each bam individually and merge the resulting peaks
    var mergedPeaks: Map<String, List<Region>>? = null
    if (pileUpInputs.size > 1) {
        for (pileUpInput in pileUpInputs) {
            val pileUps = runPileUp(listOf(pileUpInput))
            val pdfs = runSmooth(pileUps, smoothing, normalizePDF)
            val peaks = callPeaks(pdfs, threshold)
            mergedPeaks = if (mergedPeaks == null) peaks else mergePeaks(mergedPeaks, peaks)
        }
    }

    // Run peaks on all the bams piled up together, then merge them
    val pileUps = runPileUp(pileUpInputs)
    if (signalOut != null && signalOut.type == SignalOutputType.RAW) {
        createSignalFile(signalOut.path, signalOut.format, pileUps)
    }

    val pdfs = runSmooth(pileUps, smoothing, normalizePDF)
    if (signalOut != null && signalOut.type == SignalOutputType.SMOOTHED) {
        createSignalFile(signalOut.path, signalOut.format, pdfs, signalOut.signalResolution)
    }

    val peaks = callPeaks(pdfs, threshold)
    mergedPeaks = if (mergedPeaks == null) peaks else mergePeaks(mergedPeaks, peaks)

    if (peaksOut == null) return
    if (fitMode == FitMode.SKEW) {
        val skewSubPeaks = SkewFitter.fitAll(mergedPeaks, pdfs)
        writeSkewSubPeaksBed(peaksOut, skewSubPeaks)
    } else {
        val subPeaks = StandardFitter.fitAll(mergedPeaks, pdfs)
        writeStandardSubPeaksBed(peaksOut, subPeaks)
    }
     */
}

fun main(args: Array<String>) = ZPeaksCommand().main(args)
