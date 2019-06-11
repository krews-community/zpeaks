import io.SignalOutputFormat
import model.*
import mu.KotlinLogging
import org.junit.jupiter.api.*
import step.*
import step.subpeaks.*
import util.*
import java.nio.file.*

private val log = KotlinLogging.logger {}

@Disabled
class PerformanceTest {

    @Test
    fun `Run Skew Sub-Peaks on Large Peak`() {
        val pileUp = pileUp()
        val pdf = pdf(TEST_BAM_CHR, pileUp, 50.0, false)
        val peaks = callChromPeaks(pdf, 6.0)

        val maxPeak = peaks.maxBy { it.region.end - it.region.start }!!
        SkewFitter.fitChrom(TEST_BAM_CHR, listOf(maxPeak), pdf)
    }

    @Test
    fun `Run Skew Sub-Peaks on All Peaks`() {
        val testSamIn = TEST_BAM_2_PATH
            //TEST_BAM_PATH
        val signalFilename = "${testSamIn.filenameWithoutExtension()}.signal.bedGraph"
        val subPeaksFilename = "${testSamIn.filenameWithoutExtension()}.subPeaks.bed"
        val testDir = Files.createTempDirectory("zpeaks_test")
        val signalOut = testDir.resolve(signalFilename)
        val subPeaksOut = testDir.resolve(subPeaksFilename)

        run(ZPeaksRunConfig(listOf(PileUpInput(testSamIn, PileUpOptions(Strand.BOTH, PileUpAlgorithm.START))),
            signalOut = SignalOutput(signalOut, SignalOutputType.SMOOTHED, SignalOutputFormat.BED_GRAPH),
            peaksOut = null,
            subPeaksOut = subPeaksOut,
            smoothing = 50.0,
            normalizePDF = false,
            threshold = 6.0,
            fitMode = FitMode.SKEW))

        signalOut.copyToAndDelete(testSamIn.resolveSibling(signalFilename))
        subPeaksOut.copyToAndDelete(testSamIn.resolveSibling(subPeaksFilename))
    }

    @Test
    fun `Run Skew Sub-Peaks on All Peaks from Multi-File Pile-Up`() {
        val pileUpInputs =
            listOf(MULTI_BAM_1_PATH, MULTI_BAM_2_PATH, MULTI_BAM_3_PATH, MULTI_BAM_4_PATH)
                .map { PileUpInput(it, PileUpOptions(Strand.BOTH, PileUpAlgorithm.START)) }
        val subPeaksFilename = "multi_test.subPeaks.bed"
        val testDir = Files.createTempDirectory("zpeaks_test")
        val subPeaksOut = testDir.resolve(subPeaksFilename)
        run(ZPeaksRunConfig(pileUpInputs,
            signalOut = null,
            peaksOut = null,
            subPeaksOut = subPeaksOut,
            smoothing = 50.0,
            normalizePDF = false,
            threshold = 6.0,
            fitMode = FitMode.SKEW))

        subPeaksOut.copyToAndDelete(MULTI_BAM_1_PATH.resolveSibling(subPeaksFilename))
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

private fun pileUp() = runPileUp(TEST_BAM_PATH, PileUpOptions(Strand.BOTH, PileUpAlgorithm.START))
    .getValue(TEST_BAM_CHR)

private fun runManyAndAverage(sampleRange: IntRange, fitter: Fitter<*>) {
    val pileUp = pileUp()

    val errors = mutableSetOf<Double>()
    val times = mutableSetOf<Long>()
    repeat(25) {
        log.info { "Running $it" }
        val pdf = pdf(TEST_BAM_CHR, pileUp, 50.0, false, sampleRange)
        val peaks = callChromPeaks(pdf, 6.0)
        val maxPeak = peaks.maxBy { it.region.end - it.region.start }!!

        val peakRegion = maxPeak.region
        val peakValues = (peakRegion.start..peakRegion.end).map { pdf[it] }

        val startTime = System.currentTimeMillis()
        val fits = fitter.fitPeak(peakValues, peakRegion.start)
        errors += fits.map { it.error }.average()
        times += System.currentTimeMillis() - startTime
    }
    log.info { "Average error: ${errors.average()}" }
    log.info { "Average time: ${times.average()}" }
}