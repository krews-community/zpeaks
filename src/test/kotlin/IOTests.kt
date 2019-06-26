import io.*
import model.*
import org.assertj.core.api.Assertions.*
import org.junit.jupiter.api.*
import step.*
import util.*
import java.nio.file.*

class IOTests {

    private fun testPileUpFormat(fileName: String, format: SignalOutputFormat, checkContents: Boolean = true) {
        var outputPath = Files.createTempDirectory("zpeaks_test").resolve(fileName)
        val pileUp =  runPileUp(TEST_BAM_PATH, CHR_22, CHR_22_SIZE, PileUpOptions(Strand.BOTH, PileUpAlgorithm.LENGTH))
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

}