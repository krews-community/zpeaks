import model.*
import mu.KotlinLogging
import org.junit.jupiter.api.*
import org.knowm.xchart.*
import step.*
import util.*
import org.knowm.xchart.style.lines.SeriesLines
import org.knowm.xchart.style.markers.SeriesMarkers
import runner.*
import step.subpeaks.*
import java.awt.Color
import kotlin.math.min

private val log = KotlinLogging.logger {}

// Sample Ranges for Primary Single Alignment Test File
val SAMPLE_RANGE_SMALL = 10_000_000 until 15_000_000 // Single curve
val SAMPLE_RANGE_SMALL_2 = 16_773_000 until 16_776_000 // Two curves, one small
val SAMPLE_RANGE_MED = 41_000_000 until 42_000_000
val SAMPLE_RANGE_MED_2 = 44_000_000 until 46_000_000
val SAMPLE_RANGE_LARGE = 46_075_000 until 46_100_000
val SAMPLE_RANGE_LARGEST = 46_050_000 until 46_075_000

@Disabled
class Plot {

    @AfterEach
    fun stop() = Thread.sleep(Long.MAX_VALUE)

    @Test
    fun `Plot Single Alignment PDF`() = plotPdf(SingleFileZRunner(singleFileConfig()), 20_890_000 until 20_910_000)

    @Test
    fun `Plot Multi Alignment Bottom-Up PDF`() = plotPdf(BottomUpZRunner(multiFileConfig()),
        30_000_000 until 30_025_000)

    @Test
    fun `Plot Big Multi Alignment Bottom-Up PDF`() = plotPdf(BottomUpZRunner(bigMultiFileConfig()),
        30_000_000 until 30_025_000)

    @Test
    fun `Plot Big Multi Alignment Bottom-Up PDF Zoomed-In`() = plotPdf(BottomUpZRunner(bigMultiFileConfig()),
        30_010_000 until 30_017_000)

    @Test
    fun `Plot Big Multi Alignment Top-Down PDF`() = plotPdf(TopDownZRunner(bigMultiFileConfig()),
            30_000_000 until 30_025_000)

    @Test
    fun `Plot Big Multi Alignment Top-Down PDF Zoomed-In`() = plotPdf(TopDownZRunner(bigMultiFileConfig()),
        30_010_000 until 30_017_000)

    @Test
    fun `Plot Skew Sub-Peaks`() = plotSubPeaks(SingleFileZRunner(singleFileConfig()),
        SAMPLE_RANGE_LARGEST, SkewFitter)

    @Test
    fun `Plot Standard Sub-Peaks`() = plotSubPeaks(SingleFileZRunner(singleFileConfig()),
        SAMPLE_RANGE_LARGEST, StandardFitter)

    @Test
    fun `Plot Skew Sub-Peaks from Multiple Alignments`() =
        plotSubPeaks(BottomUpZRunner(multiFileConfig()),30_000_000 until 30_400_000, SkewFitter)

    @Test
    fun `Plot Skew Sub-Peaks from Multiple Alignments - Top Down`() =
        plotSubPeaks(TopDownZRunner(multiFileConfig()), 30_000_000 until 30_400_000, SkewFitter)

    @Test
    fun `Plot Skew Sub-Peaks from Many Alignments - Top Down`() =
        plotSubPeaks(TopDownZRunner(bigMultiFileConfig()),30_000_000 until 30_400_000, SkewFitter)
}

private fun simpleZRunConfig(pileUpInputs: List<PileUpInput>) = ZRunConfig(
    BamPileUpRunner(pileUpInputs), null, null, null, 50.0, 6.0)

private fun singleFileConfig() = simpleZRunConfig(listOf(PileUpInput(TEST_BAM_PATH, PileUpOptions(Strand.BOTH, PileUpAlgorithm.START))))

private fun multiFileConfig() = simpleZRunConfig(
    listOf(MULTI_BAM_1_PATH, MULTI_BAM_2_PATH, MULTI_BAM_3_PATH, MULTI_BAM_4_PATH)
        .map { PileUpInput(it, PileUpOptions(Strand.BOTH, PileUpAlgorithm.START)) }
)

private fun bigMultiFileConfig() = simpleZRunConfig(
    MANY_BAM_PATHS.map { PileUpInput(it, PileUpOptions(Strand.BOTH, PileUpAlgorithm.START)) }
)

private fun plotPdf(zRunner: ZRunner, range: IntRange) {
    val displayRange = range withNSteps 1000

    val pileUp = zRunner.pileUp(CHR_22, CHR_22_SIZE, range)
    val bpUnits = BPUnits.MBP
    val pileUpChartData = bpData(displayRange, bpUnits) { pileUp[it].toDouble() }
    val pileUpChart = xyAreaChart("Pile Up", bpUnits, pileUpChartData)

    val pdf = zRunner.pdf(pileUp)
    val pdfChartData = bpData(displayRange, bpUnits) {
        (pdf[it] - pdf.background.average) / pdf.background.stdDev
    }
    val pdfChart = xyAreaChart("PDF (in Standard Deviations from Avg)", bpUnits, pdfChartData)

    val peaks = zRunner.peaks(pdf)
    val peaksChartData = regionsData(peaks, displayRange, bpUnits,
        pdfChartData.yValues.min()!!, pdfChartData.yValues.max()!!)

    val series = pdfChart.addSeries("Peaks", peaksChartData.xValues, peaksChartData.yValues)
    series.xySeriesRenderStyle = XYSeries.XYSeriesRenderStyle.StepArea
    series.fillColor = Color(100, 200, 100, 100)
    series.marker = SeriesMarkers.NONE
    series.lineStyle = SeriesLines.NONE

    SwingWrapper(listOf(pileUpChart, pdfChart)).displayChartMatrix()
}

private fun <T : GaussianParameters> plotSubPeaks(zRunner: ZRunner, range: IntRange, fitter: Fitter<T>) {
    val pileUp = zRunner.pileUp(CHR_22, CHR_22_SIZE, range)
    val pdf = zRunner.pdf(pileUp)
    val peaks = zRunner.peaks(pdf)
    val maxPeak = peaks.maxBy { it.end - it.start }
    log.info { "Max Peak: $maxPeak" }

    val bpUnits = BPUnits.KBP
    val peakRegion = maxPeak!!
    val displayRange = peakRegion.start until peakRegion.end - 1 withNSteps 1000
    val peaksChartData = bpData(displayRange, bpUnits) { pdf[it].toDouble() }
    val chart = xyAreaChart("Peak", bpUnits, peaksChartData)

    val peakValues = (maxPeak.start..maxPeak.end).map { pdf[it] }
    val fits = fitter.fitPeak(peakValues.map { it.toDouble() }, maxPeak.start)

    for ((fitIndex, fit) in fits.withIndex()) {
        val fitDisplayRange = fit.region.start until fit.region.end - 1 withNSteps 300

        val subPeakValues = fit.subPeaks.map { subPeak ->
            val params = subPeak.gaussianParameters
            fitter.optimizer.calculateCurve(listOf(params), fit.region.end - fit.region.start, fit.region.start)
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

fun regionsData(regions: List<Region>, range: IntProgression, bpUnits: BPUnits, offValue: Double, onValue: Double): XYChartData {
    val points = mutableListOf<Pair<Double, Double>>()
    points += range.first / bpUnits.bpPerUnit.toDouble() to offValue
    points += range.last / bpUnits.bpPerUnit.toDouble() to offValue

    val regionsInRange = regions.filter { it.end >= range.first || it.start <= range.last }
    regionsInRange.forEach {
        val start = if (it.start < range.first) range.first else it.start
        points += start / bpUnits.bpPerUnit.toDouble() to onValue
        val end = if (it.end > range.last) range.last else it.end
        points += end / bpUnits.bpPerUnit.toDouble() to offValue
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