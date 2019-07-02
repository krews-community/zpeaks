package util

import java.lang.Float.floatToRawIntBits
import java.lang.Float.intBitsToFloat
import java.util.concurrent.atomic.AtomicIntegerArray

class AtomicFloatArray(val size: Int) {
    @Transient
    private val ints = AtomicIntegerArray(size)

    operator fun get(i: Int): Float {
        return intBitsToFloat(ints[i])
    }

    operator fun set(i: Int, newValue: Float) {
        ints[i] = floatToRawIntBits(newValue)
    }

    fun addAndGet(i: Int, delta: Float): Float {
        while (true) {
            val current = ints.get(i)
            val currentVal = intBitsToFloat(current)
            val nextVal = currentVal + delta
            val next = floatToRawIntBits(nextVal)
            if (ints.compareAndSet(i, current, next)) {
                return nextVal
            }
        }
    }

    fun average(): Double {
        var sum = 0.0
        var count = 0
        for (element in 0 until size) {
            sum += element
            ++count
        }
        return if (count == 0) Double.NaN else sum / count
    }

}