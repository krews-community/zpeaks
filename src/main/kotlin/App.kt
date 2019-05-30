import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.*
import io.*
import model.*
import step.*
import step.subpeaks.*
import util.*
import java.nio.file.Path


class ZPeaks : CliktCommand() {
    private val samIn: Path by option("-samIn", help="Input Sam or Bam alignment file").path().required()
    private val signalOutPath: Path? by option("-signalOut", help="Output Pile-Up file").path()
    private val signalOutType: SignalOutputType by option("-signalOutType")
        .choice(SignalOutputType.values().associateBy { it.lowerHyphenName })
        .default(SignalOutputType.SMOOTHED)
    private val signalOutFormat: SignalOutputFormat by option("-signalOutFormat")
        .choice(SignalOutputFormat.values().associateBy { it.lowerHyphenName })
        .default(SignalOutputFormat.BIG_WIG)
    private val peaksOut: Path? by option("-peaksOut", help="Output Peaks bed file").path()
    private val subPeaksOut: Path? by option("-subPeaksOut", help="Output Sub-Peaks bed file").path()
    private val strand: Strand by option("-strand", help="Strand to count during pile-up")
        .choice(Strand.values().associateBy { it.lowerHyphenName })
        .default(Strand.BOTH)
    private val forwardShift by option("-forwardShift", help="During pile-up, shift the " +
            "forward strand by this amount. Can be positive or negative. Default 0.")
        .int().default(0)
    private val reverseShift by option("-reverseShift", help="During pile-up, shift the " +
            "reverse strand by this amount. Can be positive or negative. Default 0.")
        .int().default(0)
    private val pileUpAlgorithm: PileUpAlgorithm by option("-pileUpAlgorithm",
        help="Algorithm used to select values during pile up.")
        .choice(PileUpAlgorithm.values().associateBy { it.lowerHyphenName })
        .default(PileUpAlgorithm.START)
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
        if (signalOutPath == null && peaksOut == null && subPeaksOut == null)
            throw Exception("One of the following must be set: -signalOut, -peaksOut, or -subPeaksOut")
        val signalOut =
            if (signalOutPath != null) SignalOutput(signalOutPath!!, signalOutType, signalOutFormat)
            else null
        val pileUpOptions = PileUpOptions(strand, pileUpAlgorithm, forwardShift, reverseShift)
        run(samIn, signalOut, peaksOut, subPeaksOut, pileUpOptions, smoothing, normalizePDF,
            threshold, fitMode, parallelism)
    }
}

enum class FitMode { SKEW, STANDARD }
data class SignalOutput(val path: Path, val type: SignalOutputType, val format: SignalOutputFormat)
enum class SignalOutputType { RAW, SMOOTHED }

fun run(samIn: Path, signalOut: SignalOutput?, peaksOut: Path?, subPeaksOut: Path?, pileUpOptions: PileUpOptions,
        smoothing: Double, normalizePDF: Boolean, threshold: Double, fitMode: FitMode = FitMode.SKEW,
        parallelism: Int? = null) {
    if (parallelism != null) {
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", parallelism.toString())
    }

    val pileUps = runPileUp(samIn, pileUpOptions)
    if (signalOut != null && signalOut.type == SignalOutputType.RAW) {
        createSignalFile(signalOut.path, signalOut.format, pileUps)
    }

    val pdfs = runSmooth(pileUps, smoothing, normalizePDF)
    if (signalOut != null && signalOut.type == SignalOutputType.SMOOTHED) {
        createSignalFile(signalOut.path, signalOut.format, pdfs)
    }

    val peaks = callPeaks(pdfs, threshold)
    if (peaksOut != null) {
        writePeaksBed(peaksOut, peaks)
    }

    if (subPeaksOut == null) return
    if (fitMode == FitMode.SKEW) {
        val skewSubPeaks = SkewFitter.fitAll(peaks, pdfs)
        writeSkewSubPeaksBed(subPeaksOut, skewSubPeaks)
    } else {
        val subPeaks = StandardFitter.fitAll(peaks, pdfs)
        writeStandardSubPeaksBed(subPeaksOut, subPeaks)
    }
}

fun main(args: Array<String>) = ZPeaks().main(args)
