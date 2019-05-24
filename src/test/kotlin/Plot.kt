import model.*
import org.junit.jupiter.api.*
import org.knowm.xchart.*
import step.*
import util.*
import org.knowm.xchart.style.lines.SeriesLines
import org.knowm.xchart.style.markers.SeriesMarkers
import step.subpeaks.*
import java.awt.Color
import kotlin.math.min


@Disabled
class Plot {

    @AfterEach fun stop() = Thread.sleep(Long.MAX_VALUE)

    @Test
    fun `Plot Pile Up, PDF, and Peaks`() {
        val sampleRange =
            //20_890_000 until 20_910_000
            16_774_000 until 16_776_000 // 16_774_965 until 16_774_997
        val displayRange = sampleRange withNSteps 500
        val pileUp = runPileUp(TEST_BAM_PATH, PileUpOptions(Strand.BOTH, PileUpAlgorithm.START))
            .getValue(TEST_BAM_CHR)

        val bpUnits = BPUnits.MBP
        val pileUpChartData = bpData(displayRange, bpUnits) {
            pileUp[it].toDouble()
        }
        val pileUpChart = xyAreaChart("Pile Up", bpUnits, pileUpChartData)

        val pdf = pdf(TEST_BAM_CHR, pileUp,50.0, false, sampleRange)
        val pdfChartData = bpData(displayRange, bpUnits) {
            (pdf[it] - pdf.background.average) / pdf.background.stdDev
        }
        val pdfChart = xyAreaChart("PDF (in Standard Deviations from Avg)", bpUnits, pdfChartData)

        val threshold = 6.0
        val peaks = callChromPeaks(pdf, threshold)
        val peaksChartData = regionsData(peaks.map { it.region }, displayRange, bpUnits)
        val peaksChart = xyAreaChart("Peaks Over $threshold", bpUnits, peaksChartData)

        SwingWrapper(listOf(pileUpChart, pdfChart, peaksChart)).displayChartMatrix()
    }

    @Test
    fun `Plot Sub-Peaks`() {
        val sampleRange =
            //10_000_000 until 15_000_000 // Small (Single curve)
            16_773_000 until 16_776_000 // Small (Two curves, one small)
            //41_000_000 until 42_000_000 // Medium
            //44_000_000 until 46_000_000 // Medium
            //46_075_000 until 46_100_000 // Large
            //46_050_000 until 46_075_000 // Largest
        val pileUp = runPileUp(TEST_BAM_PATH, PileUpOptions(Strand.BOTH, PileUpAlgorithm.START))
            .getValue(TEST_BAM_CHR)

        val pdf = pdf(TEST_BAM_CHR, pileUp, 50.0, false, sampleRange)
        val peaks = callChromPeaks(pdf, 6.0)
        val maxPeak = peaks.maxBy { it.region.end - it.region.start }

        val bpUnits = BPUnits.KBP
        val peakRegion = maxPeak!!.region
        val displayRange = peakRegion.start until peakRegion.end-1 withNSteps 1000
        val peaksChartData = bpData(displayRange, bpUnits) { pdf[it] }
        val chart = xyAreaChart("Peak", bpUnits, peaksChartData)

        val region = maxPeak.region
        val peakValues = (region.start..region.end).map { pdf[it] }
        val fits = fitSkew(peakValues, maxPeak.region.start)

        for ((fitIndex, fit) in fits.withIndex()) {
            val fitDisplayRange = fit.region.start until fit.region.end-1 withNSteps 300

            val subPeakValues = fit.subPeaks.map { subPeak ->
                val params = subPeak.gaussianParameters
                calculateSkewCurve(
                    doubleArrayOf(params.amplitude, params.mean, params.stdDev, params.shape),
                    fit.region.end - fit.region.start, fit.region.start
                )
            }
            val subPeaksChartData = bpData(fitDisplayRange, bpUnits) { bp ->
                subPeakValues.sumByDouble { curve ->
                    curve[bp - fit.region.start]
                } + fit.background
            }
            val subPeaksSeries = chart.addSeries("Fit-$fitIndex", subPeaksChartData.xValues, subPeaksChartData.yValues)
            subPeaksSeries.lineColor = Color.BLACK
            subPeaksSeries.marker = SeriesMarkers.NONE

            for ((subPeakIndex, subPeak) in subPeakValues.withIndex()) {
                val subPeakChartData = bpData(fitDisplayRange, bpUnits) { bp ->
                    subPeak[bp - fit.region.start] + fit.background
                }
                val subPeakSeries = chart.addSeries(
                    "Fit-$fitIndex:SubPeak-$subPeakIndex",
                    subPeakChartData.xValues,
                    subPeakChartData.yValues
                )
                subPeakSeries.marker = SeriesMarkers.NONE
                subPeakSeries.lineStyle = SeriesLines.DOT_DOT
            }
        }

        SwingWrapper(chart).displayChart()
    }

}

infix fun IntProgression.withNSteps(n: Int): IntProgression {
    val length = this.last - this.first
    return this step (this.last - this.first) / min(n, length)
}

enum class BPUnits(val bpPerUnit: Int) {
    BP(1), KBP(1000), MBP(1_000_000)
}

data class XYChartData(val xValues: List<Double>, val yValues: List<Double>)

fun bpData(range: IntProgression, bpUnits: BPUnits, getValue: (Int) -> Double): XYChartData {
    val xValues = mutableListOf<Double>()
    val yValues = mutableListOf<Double>()

    for (x in range) {
        val chunkRange = x until min(range.last, x + range.step)
        val chunkValue = chunkRange.map { getValue(it) }.max() ?: 0.0
        xValues += x / bpUnits.bpPerUnit.toDouble()
        yValues += chunkValue
    }

    return XYChartData(xValues, yValues)
}

fun regionsData(regions: List<Region>, range: IntProgression, bpUnits: BPUnits): XYChartData {
    val points = mutableListOf<Pair<Double, Double>>()
    points += range.first / bpUnits.bpPerUnit.toDouble() to 0.0
    points += range.last / bpUnits.bpPerUnit.toDouble() to 0.0

    val regionsInRange = regions.filter { it.start >= range.first && it.end <= range.last }
    regionsInRange.forEach {
        points += it.start / bpUnits.bpPerUnit.toDouble() to 1.0
        points += it.end / bpUnits.bpPerUnit.toDouble() to 0.0
    }
    points.sortBy { it.first }

    return XYChartData(points.map { it.first }, points.map { it.second })
}

fun xyAreaChart(title: String, bpUnits: BPUnits, data: XYChartData): XYChart {
    val chart = XYChartBuilder()
        .title(title)
        .width(800)
        .height(600)
        .xAxisTitle("Base Pair (in $bpUnits)")
        .build()
    val series = chart.addSeries("Data", data.xValues, data.yValues)
    series.xySeriesRenderStyle = XYSeries.XYSeriesRenderStyle.StepArea
    series.marker = SeriesMarkers.NONE
    series.lineStyle = SeriesLines.NONE
    chart.styler.isLegendVisible = false

    return chart
}