package statalign.base;

/**
 * Dynamically growing circular array with extremely efficient deque and hash operations.
 * 
 * Can be used as a dynamically growing hash (with semi-contiguous key range for
 * space-efficiency), deque (queue or stack), array, perhaps in turns.
 * Offers amortised O(1) time for each implemented operation.
 * Does not shrink.
 *
 * @author novadam
 *
 * @param <E> element type
 */
public class CircularArray<E> {

	/** Minimum capacity of the array, must be a power of two */
	private final static int MIN_CAPACITY = 4;

	/** lowest key of elements in array */
	int startKey;
	/** highest key of elements in array plus 1 */
	int endKey;

	/** data array, capacity of the <tt>CircularArray</tt> is <tt>data.length</tt> */
	E[] data;

	/**
	 * Constructs an empty <tt>CircularArray</tt> with an initial capacity of <i>MIN_CAPACITY</i>.
	 */
	@SuppressWarnings("unchecked")
	public CircularArray() {
		data = (E[]) new Object[MIN_CAPACITY];
	}
	
	/**
	 * Constructs an empty <tt>CircularArray</tt> with a given minimal initial capacity
	 * but at least <i>MIN_CAPACITY</i>.
	 *  
	 * @param initialCapacity  minimum capacity of the new array
	 */
	@SuppressWarnings("unchecked")
	public CircularArray(int initialCapacity) {
		if(initialCapacity > MIN_CAPACITY)
			initialCapacity = findCapacity(initialCapacity-1);
		else
			initialCapacity = MIN_CAPACITY;
		data = (E[]) new Object[initialCapacity];
	}

	/**
	 * Difference between the highest and lowest key of elements in the array plus 1.
	 * This is the length of the data area which is considered nonempty.
	 * If array is used solely as a deque (or at least keys are contiguous), this equals the
	 * number of elements in the array.
	 * 
	 * @return  length as defined above
	 */
	public int length() {
		return endKey-startKey;
	}

	/**
	 * Puts a new (key,element) pair into the <tt>CircularArray</tt>, growing the array if
	 * the new key range doesn't fit into the array.
	 * 
	 * This operation might destroy the property of contiguous keys.
	 * 
	 * @param key  key of the element
	 * @param element  the element to store at the given key
	 */
	public void put(int key, E element) {
		int sizem1 = data.length-1;
		final int sk = startKey;
		final int ek = endKey;
		if(sk == ek) {							// empty array
			endKey = (startKey = key) + 1;
		} else if(key < sk) {				// insertion before first key, length increases
			final int newSizem1 = ek-key-1;
			if(newSizem1 > sizem1) {			// not enough space, must grow array
				sizem1 = grow(newSizem1)-1;
			}
			startKey = key;
		} else if(key >= ek) {			// insertion after last key, length increases
			final int newSizem1 = key-sk;
			if(newSizem1 > sizem1) {		// not enough space, must grow array
				sizem1 = grow(newSizem1)-1;					// no need to add 1, see grow
			}
			endKey = key + 1;
		}
		data[key & sizem1] = element;
	}

	/**
	 * Retrieves the element at the given key in the array or <tt>null</tt> if there is no element
	 * with that key.
	 * 
	 * @param key  the given key
	 * @return  element at the key or <tt>null</tt>
	 */
	public E get(int key) {
		if(key < startKey || key >= endKey)			// out of stored key range
			return null;
		return data[key & (data.length-1)];
	}

	/**
	 * Adds new element at the end of the array (key <tt>endKey</tt> is assigned), growing the
	 * array if full, then <tt>endKey</tt> is increased.
	 * 
	 * If key range is contiguous before the call, this operation will maintain the property. 
	 * 
	 * @param elem  the element to be added
	 */
	public void push(E elem) {
		int size = data.length;
		if(endKey - startKey == size)
			grow2(size <<= 1);
		data[endKey++ & (size - 1)] = elem;
	}

	/**
	 * Retrieves the last element in the array (element with key <tt>endKey-1</tt>) and
	 * decreases <tt>endKey</tt>.
	 * 
	 * If key range is contiguous, a <tt>null</tt> value can only be returned if array is empty or
	 * if <tt>null</tt> values have been added to the array explicitly. Otherwise, <tt>null</tt>
	 * may be returned for keys with an unspecified value.
	 * 
	 * @return  the last element in the array or <tt>null</tt> if empty
	 */
	public E pop() {
		if(startKey == endKey)		// array empty
			return null;
		final int index = --endKey & (data.length - 1);
		final E result = data[index];
		data[index] = null;		// should we care ?? (only has an effect if array and deque operations used in turns)
		return result;
	}

	/**
	 * Adds new element at the beginning of the array (key <tt>startKey-1</tt> is assigned), growing
	 * the array if full, then <tt>startKey</tt> is decreased.
	 * 
	 * If key range is contiguous before the call, this operation will maintain the property. 
	 * 
	 * @param elem  the element to be added
	 */
	public void unshift(E elem) {
		int size = data.length;
		if(endKey - startKey == size)
			grow2(size <<= 1);
		data[--startKey & (size - 1)] = elem;
	}

	/**
	 * Retrieves the first element in the array (element with key <tt>startKey</tt>) and
	 * increases <tt>startKey</tt>.
	 * 
	 * If key range is contiguous, a <tt>null</tt> value can only be returned if array is empty or
	 * if <tt>null</tt> values have been added to the array explicitly. Otherwise, <tt>null</tt>
	 * may be returned for keys with an unspecified value.

	 * @return  the first element in the array or <tt>null</tt> if empty
	 */
	public E shift() {
		if(startKey == endKey)
			return null;
		final int index = startKey++ & (data.length - 1);
		final E result = data[index];
		data[index] = null;
		return result;
	}

	/**
	 * Copies contents of this CircularArray to a regular array. It must be at least of
	 * length() size.
	 *
	 * @return array
	 */
	public E[] toArray(E[] newData) {
		E[] oldData = data;
		int size = oldData.length;
		int sk = startKey;
		int n = endKey-sk;
		sk &= size-1;
		copy(oldData, newData, sk, 0, size, n);
		return newData;
	}

	/**
	 * Finds minimal capacity <tt>2^n >= size+1</tt>.
	 * 
	 * @param size  lower bound for capacity
	 * @return target capacity that is a power of 2 and above <tt>size</tt>
	 */
	private int findCapacity(int size) {
		// Find the best power of two to hold elements.
		// Code borrowed (= stolen) from java.util.ArrayDeque
		size |= (size >>>  1);
		size |= (size >>>  2);
		size |= (size >>>  4);
		size |= (size >>>  8);
		size |= (size >>> 16);
		size++;

		if (size < 0)   // Too many elements, must back off
			size >>>= 1;		// Good luck allocating 2^30 elements

		return size;
	}
	
	/**
	 * Grows array. New size will be <tt>2^n >= size+1</tt>.
	 * 
	 * @param size  lower bound for new size, must be above current size (not checked!)
	 * @return  new size
	 */
	private int grow(int size) {
		final int newSize = findCapacity(size);
		grow2(newSize);
		return newSize;
	}

	/**
	 * Grows array.
	 * 
	 * @param size  new size, must be a power of two
	 */
	@SuppressWarnings("unchecked")
	private void grow2(int newSize) {
		E[] newData = (E []) new Object[newSize];
		E[] oldData = data;
		int size = oldData.length;
		int sk = startKey;
		int n = endKey-sk;
		int nsk = sk & (newSize-1);
		sk &= size-1;
		int rem1 = newSize-nsk;
		int rem2 = n-rem1;
		if(rem2 > 0) {
			System.arraycopy(oldData, sk, newData, nsk, rem1);
			System.arraycopy(oldData, (sk+rem1) & (size-1), newData, 0, rem2);
		} else {
			copy(oldData, newData, sk, nsk, size, n);
		}
		data = newData;
	}	

	private void copy(E[] oldData, E[] newData, int sk, int nsk, int size, int n) {
		int rem1 = size-sk;
		int rem2 = n-rem1;
		if(rem2 > 0) {
			System.arraycopy(oldData, sk, newData, nsk, rem1);
			System.arraycopy(oldData, 0, newData, nsk+rem1, rem2);
		} else {
			System.arraycopy(oldData, sk, newData, nsk, n);
		}
	}


}
