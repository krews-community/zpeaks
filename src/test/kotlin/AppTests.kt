import model.*
import mu.KotlinLogging
import org.assertj.core.api.Assertions.*
import org.junit.jupiter.api.*
import step.*
import util.*
import java.nio.file.*

private val log = KotlinLogging.logger {}

class AppTests {

    @Test
    fun `Test App Run`() {
        val pileUps = runPileUp(TEST_BAM_PATH, Strand.BOTH, false, PileUpAlgorithm.START)

        val peaksFilename = "ENCFF375IJW.chr22.peaks.bed"
        val subPeaksFilename = "ENCFF375IJW.chr22.subPeaks.bed"
        val testDir = Files.createTempDirectory("zpeaks_test")
        val peaksOut = testDir.resolve(peaksFilename)
        val subPeaksOut = testDir.resolve(subPeaksFilename)
        runPeaks(pileUps, 50.0, true,6.0, peaksOut, subPeaksOut, (10_000_000..15_000_000))
        Files.copy(peaksOut, TEST_BAM_PATH.resolveSibling(peaksFilename), StandardCopyOption.REPLACE_EXISTING)
        Files.copy(subPeaksOut, TEST_BAM_PATH.resolveSibling(subPeaksFilename), StandardCopyOption.REPLACE_EXISTING)

        // Check the Peaks output file.
        val (peaksLineCount, peaksLineSample) = sampleOutputFile(peaksOut)
        val peaksLineValues = peaksLineSample!!.split("\t")
        assertThat(peaksLineValues[0]).isEqualTo("chr22")
        assertThatCode { peaksLineValues[1].toInt() }.doesNotThrowAnyException()
        assertThatCode { peaksLineValues[2].toInt() }.doesNotThrowAnyException()
        assertThat(peaksLineValues[3]).isNotNull()
        assertThatCode { peaksLineValues[4].toDouble() }.doesNotThrowAnyException()
        assertThat(peaksLineCount).isGreaterThan(0)

        // Check the Sub-Peaks output file
        val (subPeaksLineCount, subPeaksLineSample) = sampleOutputFile(subPeaksOut)
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