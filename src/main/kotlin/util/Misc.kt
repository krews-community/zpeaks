package util

import com.google.common.base.CaseFormat
import mu.KotlinLogging
import org.apache.commons.math3.util.FastMath
import java.io.DataInputStream
import java.nio.Buffer
import java.nio.CharBuffer
import java.nio.file.Path
import java.util.concurrent.*
import java.util.concurrent.atomic.AtomicInteger

private val log = KotlinLogging.logger {}

const val SQRT2PI = 2.506628275

fun logProgress(name: String, total: Int, run: (tracker: AtomicInteger) -> Unit) {
    val progressTracker = AtomicInteger(0)
    val logExecutor = Executors.newSingleThreadScheduledExecutor()
    val logJob = if (total > 0) {
        logExecutor.scheduleAtFixedRate({
            val complete = progressTracker.get()
            val percentage = "%.2f".format(complete.toDouble() / total * 100)
            log.info { "$name - $percentage% complete." }
        }, 5, 5, TimeUnit.SECONDS)
    } else null
    run(progressTracker)
    logJob?.cancel(false)
    logExecutor.shutdown()
}

fun <T> runParallel (name: String, valueLabel: String, values: List<T>, runForValue: (value: T) -> Unit) {
    val commonPool = ForkJoinPool.commonPool()
    log.info { "$name initiated for ${values.size} $valueLabel" }

    logProgress(name, values.size) { tracker ->
        // Run each search result into a "Callable" thread object to write the results.
        val tasks = values.map { value ->
            Callable<Unit> {
                try {
                    runForValue(value)
                } catch (e: Exception) {
                    log.error(e) { "Error during $name" }
                }
                tracker.incrementAndGet()
            }
        }

        // Run all the tasks and block until complete
        commonPool.invokeAll(tasks)
    }

    // Shutdown ExecutorService to prevent application termination issues
    log.info { "$name runs complete!" }
}

/**
 * Extension function to get enum name in lower-hyphen-format
 */
val Enum<*>.lowerHyphenName: String get() = CaseFormat.UPPER_UNDERSCORE.to(CaseFormat.LOWER_HYPHEN, this.name)

fun Double.pow(x: Int) = FastMath.pow(this, x)

val IntRange.length get() = this.endInclusive - this.start + 1

const val BIGWIG_MAGIC = 0x26FC_8F88

fun isBigWig(path: Path): Boolean {
    val magic = DataInputStream(path.toFile().inputStream()).readInt()
    return magic == BIGWIG_MAGIC || magic == Integer.reverseBytes(BIGWIG_MAGIC)
}

fun isBedGraph(path: Path): Boolean {
    // If this is a binary file, we don't want to read too much
    val head = CharBuffer.allocate(1000)
    path.toFile().bufferedReader().read(head)
    (head as Buffer).position(0)
    head.lineSequence().forEach {
        if (it.startsWith("track")) {
            return@forEach
        }
        val split = it.split("\\s".toRegex())
        if (split.size != 4) {
            return false
        }
        return split[1].toIntOrNull() != null && split[2].toIntOrNull() != null && split[3].toFloatOrNull() != null
    }
    return false
}