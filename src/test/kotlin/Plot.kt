import model.Region
import model.Strand
import org.junit.jupiter.api.*
import org.knowm.xchart.*
import step.*
import util.getResourcePath
import org.knowm.xchart.internal.chartpart.Chart
import org.knowm.xchart.style.lines.SeriesLines
import org.knowm.xchart.style.markers.SeriesMarkers
import kotlin.math.min

@Disabled
class Plot {

    @AfterEach fun stop() = Thread.sleep(Long.MAX_VALUE)

    @Test
    fun `Pile Up and PDF`() {
        val testBamPath = getResourcePath("ENCFF375IJW.chr22.bam")
        val pileUps = pileUpSam(testBamPath, Strand.BOTH, false, PileUpAlgorithm.START)
        val chr22PileUp = pileUps.getValue("chr22")

        val displayRange = 20_890_000 until 20_910_000 step 100
        val filteredPileUp = chr22PileUp.values
            .filter { it.key >= displayRange.first && it.key <= displayRange.last }

        val pileUpChart = bpDataChart("Pile Up", displayRange, BPUnits.MBP) {
            filteredPileUp.getOrDefault(it, 0).toDouble()
        }

        val pdf = pdf("chr22", filteredPileUp, chr22PileUp.sum, chr22PileUp.chromosomeLength,
            50.0, false)
        val pdfChart = bpDataChart("PDF", displayRange, BPUnits.MBP) {
            (pdf.values.getOrDefault(it, 0.0) - pdf.background.average) / pdf.background.stdDev
        }

        val threshold = 6.0
        val peaks = callPeaks(pdf, threshold)
        val peaksChart = regionsChart("Peaks Over $threshold", peaks, displayRange, BPUnits.MBP)

        SwingWrapper(listOf(pileUpChart, pdfChart, peaksChart)).displayChartMatrix()
    }

}

enum class BPUnits(val bpPerUnit: Int) {
    BP(1), KBP(1000), MBP(1_000_000)
}

fun bpDataChart(title: String, range: IntProgression, bpUnits: BPUnits, getValue: (Int) -> Double): Chart<*, *> {
    val xValues = mutableListOf<Double>()
    val yValues = mutableListOf<Double>()

    for (x in range) {
        val chunkRange = x until min(range.last, x + range.step)
        val chunkValue = chunkRange.map { getValue(it) }.max() ?: 0.0
        xValues += x / bpUnits.bpPerUnit.toDouble()
        yValues += chunkValue
    }

    return xyAreaChart(title, bpUnits, xValues, yValues)
}

fun regionsChart(title: String, regions: List<Region>, range: IntProgression, bpUnits: BPUnits): Chart<*, *> {
    val points = mutableListOf<Pair<Double, Double>>()
    points += range.first / bpUnits.bpPerUnit.toDouble() to 0.0
    points += range.last / bpUnits.bpPerUnit.toDouble() to 0.0

    val regionsInRange = regions.filter { it.start >= range.first && it.end <= range.last }
    regionsInRange.forEach {
        points += it.start / bpUnits.bpPerUnit.toDouble() to 1.0
        points += it.end / bpUnits.bpPerUnit.toDouble() to 0.0
    }
    points.sortBy { it.first }

    return xyAreaChart(title, bpUnits, points.map { it.first }, points.map { it.second })
}

fun xyAreaChart(title: String, bpUnits: BPUnits, xValues: List<Number>, yValues: List<Number>): Chart<*, *> {
    val chart = XYChartBuilder()
        .title(title)
        .width(800)
        .height(600)
        .xAxisTitle("Base Pair (in $bpUnits)")
        .build()
    val series = chart.addSeries("Data", xValues, yValues)
    series.xySeriesRenderStyle = XYSeries.XYSeriesRenderStyle.StepArea
    series.marker = SeriesMarkers.NONE
    series.lineStyle = SeriesLines.NONE
    chart.styler.isLegendVisible = false

    return chart
}