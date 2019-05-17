package util

import com.google.common.base.CaseFormat
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath
import java.util.concurrent.Executors
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicInteger

private val log = KotlinLogging.logger {}

const val SQRT2PI = 2.506628275

fun logProgress(name: String, total: Int, run: (tracker: AtomicInteger) -> Unit) {
    val progressTracker = AtomicInteger(0)
    val logExecutor = Executors.newSingleThreadScheduledExecutor()
    val logJob = logExecutor.scheduleAtFixedRate({
        val complete = progressTracker.get()
        val percentage = "%.2f".format(complete.toDouble() / total * 100)
        log.info { "$name - $percentage% complete." }
    }, 5, 5, TimeUnit.SECONDS)
    run(progressTracker)
    logJob?.cancel(false)
    logExecutor.shutdown()
}

/**
 * Extension function to get enum name in lower-hyphen-format
 */
val Enum<*>.lowerHyphenName: String get() = CaseFormat.UPPER_UNDERSCORE.to(CaseFormat.LOWER_HYPHEN, this.name)

fun Double.pow(x: Int) = FastMath.pow(this, x)