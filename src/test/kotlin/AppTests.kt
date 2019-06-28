import mu.KotlinLogging
import org.assertj.core.api.Assertions.*
import org.junit.jupiter.api.*
import java.nio.file.*
import step.*
import step.subpeaks.*
import util.*
import io.*
import model.*
import runner.SingleFileZRunner
import runner.ZRunConfig

private val log = KotlinLogging.logger {}

class AppTests {

    @Test
    fun `Test App Run`() {
        val peaksFilename = "ENCFF375IJW.chr22.peaks.bed"
        val testDir = Files.createTempDirectory("zpeaks_test")
        var peaksOut = testDir.resolve(peaksFilename)

        val runConfig = ZRunConfig(listOf(PileUpInput(TEST_BAM_PATH, PileUpOptions(Strand.BOTH, PileUpAlgorithm.START))))
        val zRunner = SingleFileZRunner(runConfig)
        val pileUp = zRunner.pileUp(CHR_22, CHR_22_SIZE, 10_000_000 until 15_000_000)

        val pdf = zRunner.pdf(pileUp)
        val peaks = zRunner.peaks(pdf)
        val subPeaks = SkewFitter.fit(CHR_22, peaks, pdf)
        writeSkewSubPeaksBed(peaksOut, CHR_22, subPeaks)

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