import io.*
import model.*
import org.assertj.core.api.Assertions.*
import org.junit.jupiter.api.*
import runner.SingleFileZRunner
import runner.ZRunConfig
import step.*
import util.*
import java.nio.file.*

class IOTests {

    private fun testPileUpFormat(fileName: String, format: SignalOutputFormat, checkContents: Boolean = true) {
        var outputPath = Files.createTempDirectory("zpeaks_test").resolve(fileName)
        prepBams(listOf(TEST_BAM_PATH))
        val pileUp =  runPileUp(TEST_BAM_PATH, CHR_22, CHR_22_SIZE,
            0 until CHR_22_SIZE, PileUpOptions(Strand.BOTH, PileUpAlgorithm.LENGTH))
        createSignalFile(outputPath, format, CHR_22, pileUp)

        outputPath = outputPath.copyToAndDelete(TEST_BAM_PATH.resolveSibling(fileName))
        if (checkContents) assertThat(outputPath).hasSameContentAs(getResourcePath(fileName))
    }

    @Test
    fun `test pile-up output in wig format`() {
        testPileUpFormat("ENCFF375IJW.chr22.pileup.wig", SignalOutputFormat.WIG)
    }

    @Test
    fun `test pile-up output in bed-graph format`() {
        testPileUpFormat("ENCFF375IJW.chr22.pileup.bedGraph", SignalOutputFormat.BED_GRAPH)
    }

    @Test
    fun `Test raw pile-up bed-graph from multiple chromosomes`() {
        val testSamIn = TEST_BAM_2_PATH
        val signalFilename = "${testSamIn.filenameWithoutExtension()}.signal.bedGraph"
        val testDir = Files.createTempDirectory("zpeaks_test")
        val signalOut = testDir.resolve(signalFilename)

        val runConfig = ZRunConfig(
            pileUpInputs = listOf(PileUpInput(testSamIn, PileUpOptions(Strand.BOTH, PileUpAlgorithm.START))),
            chrFilter = listOf("chr1", "chr2", "chr3").map { it to null }.toMap(),
            signalOut = SignalOutput(signalOut, SignalOutputType.RAW, SignalOutputFormat.BED_GRAPH),
            fitMode = FitMode.SKEW
        )
        SingleFileZRunner(runConfig).run()

        val chromosomesOut = Files.readAllLines(signalOut)
            .map { if (it.startsWith("track")) null else it.split("\t")[0] }
            .filter { it != null }
            .distinct()
        assertThat(chromosomesOut).hasSize(3)
    }

}