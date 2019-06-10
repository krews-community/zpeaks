import model.Strand
import org.assertj.core.api.Assertions.*
import org.junit.jupiter.api.*
import step.*
import util.*

@Disabled
class PileUpTests {

    @Test
    fun `Test pile-up on multiple alignment files`() {
        val pileUpOptions = PileUpOptions(Strand.BOTH, PileUpAlgorithm.START)

        val pileUpsA = runPileUp(MULTI_BAM_1_PATH, pileUpOptions)
        assertThat(pileUpsA).containsKey("chr22")
        val maxBetween0And10mA = pileUpsA.getValue("chr22").maxBetween(0, 10_000_000)
        assertThat(maxBetween0And10mA).isEqualTo(0)
        val maxBetween30And32mA = pileUpsA.getValue("chr22").maxBetween(30_000_000, 32_000_000)
        assertThat(maxBetween30And32mA).isGreaterThan(0)

        val multiPileUpInputs =
            listOf(MULTI_BAM_1_PATH, MULTI_BAM_2_PATH, MULTI_BAM_3_PATH, MULTI_BAM_4_PATH)
            .map { PileUpInput(it, pileUpOptions) }
        val pileUpsB = runPileUp(multiPileUpInputs)
        assertThat(pileUpsB).containsKey("chr22")
        val maxBetween0And10mB = pileUpsB.getValue("chr22").maxBetween(0, 10_000_000)
        assertThat(maxBetween0And10mB).isEqualTo(0)
        val maxBetween30And32mB = pileUpsB.getValue("chr22").maxBetween(30_000_000, 32_000_000)
        assertThat(maxBetween30And32mB).isGreaterThan(0)

        assertThat(maxBetween30And32mB).isGreaterThan(maxBetween30And32mA)
    }

}

fun PileUp.maxBetween(from: Int, to: Int): Int {
    var max = 0
    for (i in from until to) {
        if (this[i] > max) max = this[i]
    }
    return max
}