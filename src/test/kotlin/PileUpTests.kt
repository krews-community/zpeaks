import model.*
import org.assertj.core.api.Assertions.*
import org.junit.jupiter.api.*
import runner.BottomUpZRunner
import runner.ZRunConfig
import step.*
import util.*

@Disabled
class PileUpTests {

    @Test
    fun `Test pile-up on multiple alignment files`() {
        val pileUpOptions = PileUpOptions(Strand.BOTH, PileUpAlgorithm.START)

        val pileUpsA = runPileUp(MULTI_BAM_1_PATH, CHR_22, CHR_22_SIZE, pileUpOptions)
        assertThat(pileUpsA.chr).isEqualTo(CHR_22)
        val maxBetween0And10mA = pileUpsA.maxBetween(0, 10_000_000)
        assertThat(maxBetween0And10mA).isEqualTo(0)
        val maxBetween30And32mA = pileUpsA.maxBetween(30_000_000, 32_000_000)
        assertThat(maxBetween30And32mA).isGreaterThan(0.0)

        val multiPileUpInputs =
            listOf(MULTI_BAM_1_PATH, MULTI_BAM_2_PATH, MULTI_BAM_3_PATH, MULTI_BAM_4_PATH)
            .map { PileUpInput(it, pileUpOptions) }
        val multiRunner = BottomUpZRunner(ZRunConfig(multiPileUpInputs))
        val pileUpsB = multiRunner.pileUp(CHR_22, CHR_22_SIZE)
        assertThat(pileUpsB.chr).isEqualTo(CHR_22)
        val maxBetween0And10mB = pileUpsB.maxBetween(0, 10_000_000)
        assertThat(maxBetween0And10mB).isEqualTo(0)
        val maxBetween30And32mB = pileUpsB.maxBetween(30_000_000, 32_000_000)
        assertThat(maxBetween30And32mB).isGreaterThan(0.0)

        assertThat(maxBetween30And32mB).isGreaterThan(maxBetween30And32mA)
    }

}

fun PileUp.maxBetween(from: Int, to: Int): Double {
    var max = 0.0
    for (i in from until to) {
        if (this[i] > max) max = this[i]
    }
    return max
}