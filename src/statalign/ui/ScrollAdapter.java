package statalign.ui;

import java.awt.Window;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;

import javax.swing.JScrollBar;

public class ScrollAdapter extends KeyAdapter {
	
	private JScrollBar scroll;
	private Window window;
	
	private int scrollMultiplier = 1;

	ScrollAdapter(JScrollBar scroll) {
		this.scroll = scroll;
	}
	
	ScrollAdapter(JScrollBar scroll, int scrollMultiplier) {
		this.scroll = scroll;
		this.scrollMultiplier = scrollMultiplier;
	}
	
	ScrollAdapter(JScrollBar scroll, int scrollMultiplier, Window windowDispose) {
		this.scroll = scroll;
		this.scrollMultiplier = scrollMultiplier;
		this.window = windowDispose;
	}
	
	public int getScrollMultiplier() {
		return scrollMultiplier;
	}
	
	public void setScrollMultiplier(int scrollMultiplier) {
		this.scrollMultiplier = scrollMultiplier;
	}
	
	public void setWindowDispose(Window window) {
		this.window = window;
	}

	@Override
	public void keyPressed(KeyEvent e) {
//		System.out.println("before: "+scroll.getValue());
		if(e.getKeyCode() == KeyEvent.VK_DOWN) {
			int value = Math.min(scroll.getMaximum(),
					scroll.getValue()+scroll.getBlockIncrement()*scrollMultiplier);
//			System.out.println(value);
			scroll.setValue(value);
//			System.out.println(scroll.getValue());
		} else if(e.getKeyCode() == KeyEvent.VK_UP) {
			int value = Math.max(scroll.getMinimum(),
					scroll.getValue()-scroll.getBlockIncrement()*scrollMultiplier);
//			System.out.println(value);
			scroll.setValue(value);
//			System.out.println(scroll.getValue());
		} else if(e.getKeyCode() == KeyEvent.VK_HOME) {
			int value = scroll.getMinimum();
//			System.out.println(value);
			scroll.setValue(value);
//			System.out.println(scroll.getValue());
		} else if(e.getKeyCode() == KeyEvent.VK_END) {
			int value = scroll.getMaximum();
//			System.out.println(value);
			scroll.setValue(value);
//			System.out.println(scroll.getValue());
		} else if(e.getKeyCode() == KeyEvent.VK_ESCAPE && window != null) {
			window.dispose();
		}
	}

}
