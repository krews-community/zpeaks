import mu.KotlinLogging
import org.assertj.core.api.Assertions.*
import org.junit.jupiter.api.*
import java.nio.file.*
import step.*
import step.subpeaks.*
import util.*
import io.*
import model.*

private val log = KotlinLogging.logger {}

class AppTests {

    @Test
    fun `Test App Run`() {
        val peaksFilename = "ENCFF375IJW.chr22.peaks.bed"
        val testDir = Files.createTempDirectory("zpeaks_test")
        var peaksOut = testDir.resolve(peaksFilename)

        val pileUp = runPileUp(TEST_BAM_PATH, PileUpOptions(Strand.BOTH, PileUpAlgorithm.START))
            .getValue(TEST_BAM_CHR)

        val pdf = pdf(pileUp, 50.0, false, 10_000_000 .. 15_000_000)
        val peaks = callChromPeaks(pdf, 6.0)
        val subPeaks = SkewFitter.fitChrom(TEST_BAM_CHR, peaks, pdf)
        writeSkewSubPeaksBed(peaksOut, mapOf(TEST_BAM_CHR to subPeaks))

        peaksOut = peaksOut.copyToAndDelete(TEST_BAM_PATH.resolveSibling(peaksFilename))

        // Check the Peaks output file
        val (subPeaksLineCount, subPeaksLineSample) = sampleOutputFile(peaksOut)
        val subPeaksLineValues = subPeaksLineSample!!.split("\t")
        assertThat(subPeaksLineValues[0]).isEqualTo("chr22")
        assertThatCode { subPeaksLineValues[1].toInt() }.doesNotThrowAnyException()
        assertThatCode { subPeaksLineValues[2].toInt() }.doesNotThrowAnyException()
        assertThat(subPeaksLineValues[3]).isNotNull()
        val splitSubPeaksValue4 = subPeaksLineValues[4].split("#")
        assertThatCode { splitSubPeaksValue4[0].toDouble() }.doesNotThrowAnyException()
        assertThatCode { splitSubPeaksValue4[1].toDouble() }.doesNotThrowAnyException()
        assertThat(subPeaksLineCount).isGreaterThan(0)
    }

}

data class OutputFileSample(val lineCount: Int, val lineSample: String?)

private fun sampleOutputFile(path: Path): OutputFileSample {
    var lineCount = 0
    var lineSample: String? = null
    Files.newBufferedReader(path).use { reader ->
        reader.forEachLine { line ->
            if (0 == lineCount) lineSample = line
            lineCount++
        }
    }
    return OutputFileSample(lineCount, lineSample)
}