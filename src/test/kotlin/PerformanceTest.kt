import io.writeSkewSubPeaksBed
import model.*
import mu.KotlinLogging
import org.junit.jupiter.api.*
import step.*
import step.subpeaks.*
import util.TEST_BAM_PATH
import java.nio.file.*

private val log = KotlinLogging.logger {}

@Tag("manual")
class PerformanceTest {

    @Test
    fun `Run Skew Sub-Peaks on Large Peak`() {
        val pileUps = pileUpSam(TEST_BAM_PATH, Strand.BOTH, false, PileUpAlgorithm.START)
        val chr22PileUp = pileUps.getValue("chr22")

        val pdf = pdf("chr22", chr22PileUp, 50.0, false)
        val peaks = callPeaks(pdf, 6.0)

        val maxPeak = peaks.maxBy { it.region.end - it.region.start }
        runSkewSubPeaks(pdf, listOf(maxPeak!!), 1)
    }

    @Test
    fun `Run Skew Sub-Peaks on All Peaks`() {
        val pileUps = pileUpSam(TEST_BAM_PATH, Strand.BOTH, false, PileUpAlgorithm.START)
        val chr22PileUp = pileUps.getValue("chr22")

        val pdf = pdf("chr22", chr22PileUp, 50.0, false)
        val peaks = callPeaks(pdf, 6.0)
        val subPeaks = runSkewSubPeaks(pdf, peaks, 6)

        val subPeaksFilename = "ENCFF375IJW.chr22.subPeaks.bed"
        val testDir = Files.createTempDirectory("zpeaks_test")
        val subPeaksOut = testDir.resolve(subPeaksFilename)
        writeSkewSubPeaksBed(subPeaksOut, "chr22", subPeaks)
        Files.copy(subPeaksOut, TEST_BAM_PATH.resolveSibling(subPeaksFilename), StandardCopyOption.REPLACE_EXISTING)
    }

    @Test
    fun `Run Skew Sub-Peaks many times for one peak and take the average score`() {
        val sampleRange =
            //10_000_000 until 15_000_000 // Small - Single curve
            41_000_000 until 42_000_000 // Medium 1
            //44_000_000 until 46_000_000 // Medium 2
            //46_000_000 until 47_000_000 // Largest

        val pileUp = pileUpSam(TEST_BAM_PATH, Strand.BOTH, false, PileUpAlgorithm.START)
            .getValue(TEST_BAM_CHR)

        val errors = mutableSetOf<Double>()
        val times = mutableSetOf<Long>()
        repeat(25) {
            val pdf = pdf(TEST_BAM_CHR, pileUp, 50.0, false, sampleRange)
            val peaks = callPeaks(pdf, 6.0)
            val maxPeak = peaks.maxBy { it.region.end - it.region.start }

            val peakRegion = maxPeak!!.region
            val peakValues = (peakRegion.start..peakRegion.end).map { pdf[it] }

            val startTime = System.currentTimeMillis()
            val fits = fitSkew(peakValues, peakRegion.start)
            errors += fits.map { it.error }.average()
            times += System.currentTimeMillis() - startTime
        }
        log.info { "Average error: ${errors.average()}" }
        log.info { "Average time: ${times.average()}" }
    }

}