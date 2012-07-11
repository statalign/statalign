package statalign.base.thread;

/**
 * Ancestor for classes that can define potential "stopping points" in (some of)
 * their methods when called from within a StoppableThread. A StoppableThread can
 * be safely stopped by another thread when execution reaches such a point.
 * 
 * "Pausing points" can also be defined to allow the StoppableThread to be paused at
 * certain points.
 * 
 * @author novak
 */
public abstract class Stoppable {
	StoppableThread sthread;
	
	/**
	 * Default constructor. Finds out if the object has been created from within a
	 * StoppableThread.
	 */
	public Stoppable() {
		Thread current = Thread.currentThread();
		if(current instanceof StoppableThread) {
			sthread = (StoppableThread) current;
		}
	}
	
	/**
	 * Should be called from methods of the descendant classes at points where execution
	 * is allowed to be suspended by an external thread (by calling the StoppableThread's
	 * suspendSoft() method)
	 * 
	 * If the descendant object has not been created from within a StoppableThread, calling
	 * this method has no effect.
	 */
	public void pausable() {
		if(sthread != null)
			sthread.pausable();
	}
	/**
	 * Should be called from methods of the descendant classes at points where execution
	 * is allowed to be stopped by an external thread (by calling the StoppableThread's
	 * stopSoft() method)
	 * If this happens, a StoppedException is thrown that can be caught at any point within
	 * the StoppableThread to do things up before exiting.
	 * 
	 * If the descendant object has not been created from within a StoppableThread, calling
	 * this method has no effect.

	 * @throws StoppedException When stopping the execution of the StoppableThread has been
	 *  requested via stopSoft()
	 */
	public void stoppable() throws StoppedException {
		if(sthread != null)
			sthread.stoppable();
	}

}