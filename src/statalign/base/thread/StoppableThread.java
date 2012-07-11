package statalign.base.thread;

/**
 * Objects of this class can be used in place of a Thread object if it is intended that
 * they can be stopped or suspended "softly" by an external thread. This means that if an
 * external thread requests termination of this thread (by calling stopSoft()),
 * execution will continue up to the first "stopping point", i.e. the first call of
 * stoppable() from within this thread.
 * At this point, a StoppedException will be thrown, and the StoppableThread will be given
 * the chance to catch it at an arbitrary point and do up before termination.
 * 
 * In every other aspect, a StoppableThread is like a Thread, e.g. you should override
 * run() and place your "softly stoppable" code there. But you must declare your run()
 * method <code>synchronized</code> for most of the above to work properly.
 * 
 * The easiest way to define "stopping points" and "pausing points" in your own classes
 * is to use Stoppable as the ancestor and call its stoppable() and pausable() methods from
 * your methods.
 * 
 * @author novak
 *
 */
public class StoppableThread extends Thread {
	volatile boolean pausing;
	volatile boolean stopping;

	/**
	 * Should be called from within the thread to define a "pausing point". If an external
	 * thread has previously requested the suspension of the thread (by calling suspendSoft()),
	 * this method will block until the external thread calls resumeSoft() or stopSoft().
	 */
	void pausable() {
		while(pausing)
			try {
				wait();
			} catch(InterruptedException e) {
			}
	}
	
	/**
	 * Should be called from within the thread to define a "stopping point". This implies
	 * a "pausing point" (see pausable()). If an external thread has previously requested the
	 * termination of the thread (by calling stopSoft()), a StoppedException will be thrown
	 * which can be caught at any point within this thread to make up before termination.
	 */
	void stoppable() throws StoppedException {
		pausable();
		if(stopping)
			throw new StoppedException();
	}
	
	/**
	 * Can be called by an external thread to request suspension of this thread. The
	 * method will not return until execution of this thread reaches a "pausing" or
	 * "stopping point", i.e. until pausable() or stoppable() is called from within this
	 * thread. At that point, this thread will be suspended until resumeSoft() or stopSoft()
	 * is called by the external thread.
	 * 
	 * Compare deprecated "hard" suspend: Thread.suspend()
	 */
	public void suspendSoft() {
		pausing = true;
		synchronized(this) {}
	}
	
	/**
	 * Calls suspendSoft() or resumeSoft(), depending on the current state of this thread.
	 */
	public void pauseToggleSoft() {
		if(pausing)
			resumeSoft();
		else
			suspendSoft();
	}
	
	/**
	 * Can be called by an external thread to request suspension of this thread. Unlike
	 * pauseSoft(), this method is non-blocking, but otherwise has the same effect.
	 */
	public void pauseNoWait() {
		pausing = true;
	}
	
	/**
	 * Resumes execution of this thread if it has been previously paused by
	 * pauseSoft(). To be called by an external thread. Non-blocking, even if pauseNoWait()
	 * was used and this thread has not yet reached a "pausing point".
	 * 
	 * Compare deprecated "hard" resume: Thread.resume()
	 */
	public synchronized void resumeSoft() {
		pausing = false;
		notify();
	}
	
	/**
	 * Can be called by an external thread to request termination of this thread. This
	 * method will not return until execution of this thread reaches a 
	 * "stopping point", i.e. until stoppable() is called from within this
	 * thread.
	 * Can also be used if this thread has been suspended by suspendSoft().
	 * 
	 * Compare deprecated "hard" stop: Thread.stop()
	 */
	public void stopSoft() {
		stopping = true;
		resumeSoft();
	}
	
	/**
	 * Similar to stopSoft() but non-blocking, and also will not wake up this thread if
	 * it is suspended. In that case, use (non-blocking) resumeSoft() first.
	 */
	public void stopNoWait() {
		pausing = false;
		stopping = true;
	}

}