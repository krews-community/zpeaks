import io.writeSkewSubPeaksBed
import model.*
import org.junit.jupiter.api.*
import step.*
import util.TEST_BAM_PATH
import java.nio.file.*

@Disabled
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

}