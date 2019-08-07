import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.*
import io.*
import model.*
import runner.ZRunConfig
import step.*
import util.*
import runner.*
import java.lang.IllegalStateException
import java.nio.file.Path


class ZPeaksCommand(val run: (RunType, ZRunConfig, Boolean) -> Unit = ::run): CliktCommand() {
    private val runType: RunType by option("-runType",
        help="Type of run. Relates to the way data from bam files are aggregated.")
        .choice(RunType.values().associateBy { it.lowerHyphenName })
        .default(RunType.BOTTOM_UP)
    private val bamsIn: List<Path> by option("-bamIn", help="Input Bam alignment file. " +
            "Cannot be used with -pileupIn.")
        .path(exists = true)
        .multiple()
    private val pileupsIn: List<Path> by option("-pileupIn", help="Input bedGraph pileup. " +
            "Cannot be used with -bamIn.")
        .path(exists = true)
        .multiple()
    private val chrFilterPath: Path? by option("-chrFilter", help="Chromosome Filter file. An optional new-line " +
            "delimited file containing the chromosome that we will analyze. Others will be ignored.")
        .path(exists = true)
    private val signalOutPath: Path? by option("-signalOut", help="Output Pile-Up file").path()
    private val signalOutType: SignalOutputType by option("-signalOutType")
        .choice(SignalOutputType.values().associateBy { it.lowerHyphenName })
        .default(SignalOutputType.SMOOTHED)
    private val signalOutFormat: SignalOutputFormat by option("-signalOutFormat")
        .choice(SignalOutputFormat.values().associateBy { it.lowerHyphenName })
        .default(SignalOutputFormat.BED_GRAPH)
    private val peaksOut: Path? by option("-peaksOut", help="Output Peaks bed file").path()
    private val strands: List<Strand> by option("-strand", help="Strand to count during pile-up. " +
            "Ignored if not using -bamIn.")
        .choice(Strand.values().associateBy { it.lowerHyphenName })
        .multiple(listOf(Strand.BOTH))
    private val forwardShifts by option("-forwardShift", help="During pile-up, shift the " +
            "forward strand by this amount. Ignored if not using -bamIn. Can be positive or negative. Default 0.")
        .int()
        .multiple(listOf(0))
    private val reverseShifts by option("-reverseShift", help="During pile-up, shift the " +
            "reverse strand by this amount. Ignored if using -bamIn. Can be positive or negative. Default 0.")
        .int()
        .multiple(listOf(0))
    private val pileUpAlgorithms: List<PileUpAlgorithm> by option("-pileUpAlgorithm",
        help="Algorithm used to select values during pile up. Ignored if not using -bamIn.")
        .choice(PileUpAlgorithm.values().associateBy { it.lowerHyphenName })
        .multiple(listOf(PileUpAlgorithm.START))
    private val signalResolution: Int by option("-signalResolution",
        help="Number of decimal places to keep in outputted signal values")
        .int().default(5)
    private val smoothing: Double by option("-smoothingFactor",
        help="Smoothing Factor for calculating PDF for Pile Up data during Peaks step.").double()
        .default(50.0)
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
        if ((bamsIn.isEmpty() && pileupsIn.isEmpty()) || (!bamsIn.isEmpty() && !pileupsIn.isEmpty()))
            throw Exception("One (and only one) of the following must be set: -bamIn or -pileupIn")
        val signalOut =
            if (signalOutPath != null) SignalOutput(signalOutPath!!, signalOutType, signalOutFormat, signalResolution)
            else null

        val chrFilter = if (chrFilterPath != null) parseChrFilter(chrFilterPath!!) else null

        val pileUpRunner: PileUpRunner = if (!bamsIn.isEmpty()) {
            val pileUpInputs = mutableListOf<PileUpInput>()
            for ((samIndex, sam) in bamsIn.withIndex()) {
                val strand = strands.getOrElse(samIndex) { strands.last() }
                val pileUpAlgorithm = pileUpAlgorithms.getOrElse(samIndex) { pileUpAlgorithms.last() }
                val forwardShift = forwardShifts.getOrElse(samIndex) { forwardShifts.last() }
                val reverseShift = reverseShifts.getOrElse(samIndex) { reverseShifts.last() }
                val pileUpOptions = PileUpOptions(strand, pileUpAlgorithm, forwardShift, reverseShift)
                pileUpInputs += PileUpInput(sam, pileUpOptions)
            }
            BamPileUpRunner(pileUpInputs)
        } else if (!pileupsIn.isEmpty()) {
            when {
                pileupsIn.all { isBigWig(it) } -> throw Exception("bigWig support is not supported yet")
                pileupsIn.all { isBedGraph(it) } -> BedGraphPileUpRunner(pileupsIn)
                else -> throw Exception("All of -pileupIn must be a bedGraph")
            }
        } else {
            // Checked at start
            throw IllegalStateException("Both bamsIn and pileupsIn are empty. (This is a bug.)")
        }

        val isSingle = if (!bamsIn.isEmpty()) {
            bamsIn.size == 1
        } else if (!pileupsIn.isEmpty()) {
            pileupsIn.size == 1
        } else {
            // Checked at start
            throw IllegalStateException("Both bamsIn and pileupsIn are empty. (This is a bug.)")
        }

        val runConfig = ZRunConfig(pileUpRunner, chrFilter, signalOut, peaksOut, smoothing,
            threshold, fitMode, parallelism)
        run(runType, runConfig, isSingle)
    }
}

enum class RunType { TOP_DOWN, BOTTOM_UP }

fun run(runType: RunType, config: ZRunConfig, isSingle: Boolean) {
    val runner = when {
        isSingle -> SingleFileZRunner(config)
        runType == RunType.BOTTOM_UP -> BottomUpZRunner(config)
        else -> TopDownZRunner(config)
    }
    runner.run()
}

fun main(args: Array<String>) = ZPeaksCommand().main(args)
