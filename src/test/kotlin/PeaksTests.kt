import org.junit.jupiter.api.Test
import model.*
import org.assertj.core.api.Assertions.assertThat
import step.*

class PeaksTests {

    @Test
    fun `Test mergePeaks`() {
        val allPeaks = listOf(
            regions(2 to 4, 5 to 6, 8 to 10, 100 to 110),
            regions(1 to 3, 5 to 10, 21 to 22, 25 to 50),
            regions(30 to 40, 35 to 60, 65 to 70, 100 to 105)
        )

        val mergedPeaks = mergePeaks(allPeaks)
        val expectedPeaks = regions(1 to 4, 5 to 10, 21 to 22, 25 to 60, 65 to 70, 100 to 110)
        assertThat(mergedPeaks).isEqualTo(expectedPeaks)
    }
}

private fun regions(vararg regionPairs: Pair<Int, Int>) = regionPairs.map { Region(it.first, it.second) }