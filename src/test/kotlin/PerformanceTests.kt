import io.*
import model.*
import mu.KotlinLogging
import org.junit.jupiter.api.*
import runner.*
import step.*
import step.subpeaks.*
import util.*
import java.nio.file.*

private val log = KotlinLogging.logger {}

@Disabled
class PerformanceTest {

    @Test
    fun `Run Skew Sub-Peaks on Large Peak`() {
        val pileUp = runPileUp(TEST_BAM_PATH, CHR_22, CHR_22_SIZE,
            0 until CHR_22_SIZE, PileUpOptions(Strand.BOTH, PileUpAlgorithm.START))
        val pdf = pdf(pileUp, 50.0)
        val peaks = callPeaks(pdf, 6.0)

        val maxPeak = peaks.maxBy { it.end - it.start }!!
        SkewFitter.fit(CHR_22, listOf(maxPeak), pdf)
    }

    @Test
    fun `Run Skew Sub-Peaks on All Peaks`() {
        val testSamIn = TEST_BAM_2_PATH
        val signalFilename = "${testSamIn.filenameWithoutExtension()}.signal.bedGraph"
        val peaksFilename = "${testSamIn.filenameWithoutExtension()}.peaks.bed"
        val testDir = Files.createTempDirectory("zpeaks_test")
        val signalOut = testDir.resolve(signalFilename)
        val peaksOut = testDir.resolve(peaksFilename)

        val runConfig = ZRunConfig(
            pileUpInputs = listOf(PileUpInput(testSamIn, PileUpOptions(Strand.BOTH, PileUpAlgorithm.START))),
            signalOut = SignalOutput(signalOut, SignalOutputType.SMOOTHED, SignalOutputFormat.BED_GRAPH),
            peaksOut = peaksOut,
            fitMode = FitMode.SKEW
        )
        SingleFileZRunner(runConfig).run()

        signalOut.copyToAndDelete(testSamIn.resolveSibling(signalFilename))
        peaksOut.copyToAndDelete(testSamIn.resolveSibling(peaksFilename))
    }

    @Test
    fun `Run Bottom-Up Skew Sub-Peaks on All Peaks from Multi-File Pile-Up`() {
        val pileUpInputs =
            listOf(MULTI_BAM_1_PATH, MULTI_BAM_2_PATH, MULTI_BAM_3_PATH, MULTI_BAM_4_PATH)
                .map { PileUpInput(it, PileUpOptions(Strand.BOTH, PileUpAlgorithm.START)) }
        val peaksFilename = "multi_test.peaks.bed"
        val testDir = Files.createTempDirectory("zpeaks_test")
        val peaksOut = testDir.resolve(peaksFilename)
        val runConfig = ZRunConfig(
            pileUpInputs = pileUpInputs,
            chrFilter = mapOf(CHR_22 to (30_000_000 until 32_000_000)),
            peaksOut = peaksOut,
            fitMode = FitMode.SKEW
        )
        BottomUpZRunner(runConfig).run()

        peaksOut.copyToAndDelete(MULTI_BAM_1_PATH.resolveSibling(peaksFilename))
    }

    @Test
    fun `Run Top-Down Skew Sub-Peaks on All Peaks from Many-File Pile-Up`() {
        val pileUpInputs =
            MANY_BAM_PATHS.map { PileUpInput(it, PileUpOptions(Strand.BOTH, PileUpAlgorithm.START)) }
        val peaksFilename = "many_test.peaks.bed"
        val testDir = Files.createTempDirectory("zpeaks_test")
        val peaksOut = testDir.resolve(peaksFilename)
        val runConfig = ZRunConfig(
            pileUpInputs = pileUpInputs,
            chrFilter = mapOf(CHR_22 to (30_000_000 until 30_500_000)),
            peaksOut = peaksOut,
            fitMode = FitMode.SKEW
        )
        TopDownZRunner(runConfig).run()

        peaksOut.copyToAndDelete(MANY_BAM_PATHS[0].resolveSibling(peaksFilename))
    }

    @Test
    fun `Run Skew Sub-Peaks many times for one peak and take the average score`() {
        runManyAndAverage(SAMPLE_RANGE_LARGEST, SkewFitter)
    }

    @Test
    fun `Run Standard Sub-Peaks many times for one peak and take the average score`() {
        runManyAndAverage(SAMPLE_RANGE_LARGEST, StandardFitter)
    }

}

private fun runManyAndAverage(range: IntRange, fitter: Fitter<*>) {
    val pileUp = runPileUp(TEST_BAM_PATH, CHR_22, CHR_22_SIZE, range,
        PileUpOptions(Strand.BOTH, PileUpAlgorithm.START))

    val errors = mutableSetOf<Double>()
    val times = mutableSetOf<Long>()
    repeat(25) {
        log.info { "Running $it" }
        val pdf = pdf(pileUp, 50.0)
        val peaks = callPeaks(pdf, 6.0)
        val maxPeak = peaks.maxBy { it.end - it.start }!!

        val peakValues = (maxPeak.start..maxPeak.end).map { pdf[it] }

        val startTime = System.currentTimeMillis()
        val fits = fitter.fitPeak(peakValues, maxPeak.start)
        errors += fits.map { it.error }.average()
        times += System.currentTimeMillis() - startTime
    }
    log.info { "Average error: ${errors.average()}" }
    log.info { "Average time: ${times.average()}" }
}