package statalign.postprocess;

import java.awt.BorderLayout;
import java.awt.CardLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;

import javax.swing.AbstractButton;
import javax.swing.ButtonGroup;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JToggleButton;
import javax.swing.JToolBar;

import statalign.base.InputData;
import statalign.postprocess.gui.treeviews.TreeView;

public abstract class TreeVisualizer extends Postprocess {

    // Variables

    private JPanel mainPanel;
    private JPanel treePanel;
    private String currentPanel;
    private ButtonGroup group;
    
    protected int noOfSamples = 0;

    protected ArrayList<TreeView> treeViews;
    private HashMap<String, TreeView> treeViewsMap;

    protected ArrayList<JComponent> toolBar;

    // Functions

    protected TreeVisualizer() {
        // TODO: Plugin configuration
        screenable = true;
        hasToolBar = true;

    }
    
    public void init(TreeView[] views) {
        // Object initialization
    	//System.out.println("Initialising TreeView");
        treeViews = new ArrayList<TreeView>();
        treeViewsMap = new HashMap<String, TreeView>();
        for (TreeView view : views) {
            addTreeView(view);
        }

        // Initialize toolbar        
        toolBar = new ArrayList<JComponent>();
        createToolBar();    
    }

    @Override
    public void beforeFirstSample(InputData input) {
        super.beforeFirstSample(input);

        if(show) {
        	Enumeration<AbstractButton> buttons = group.getElements();
        	while (buttons.hasMoreElements()) {
        		buttons.nextElement().setEnabled(true);
        	}

        	// Notifies each tree view of a MCMC run.
        	for (TreeView view : treeViews) {
        		view.beforeFirstSample(input.seqs.size());
        		view.revalidate();
        	}
        }
    }

    @Override
    public JPanel getJPanel() {
        mainPanel = new JPanel(new BorderLayout());
        treePanel = new JPanel(new CardLayout());
        for (TreeView view : treeViews) {
            // Constructs a JScrollPane for each of the tree views.
            JScrollPane parent = new JScrollPane();
            parent.setViewportView(view);
            view.setParent(parent);
            treePanel.add(parent, view.getIdentifier());
        }
        mainPanel.add(treePanel, BorderLayout.CENTER);
        return mainPanel;
    }

    @Override
    public void setSampling(boolean enabled) {
        sampling = enabled;
    }

    @Override
    public ArrayList<JComponent> getToolBarItems() {
        return toolBar;
    }

    public void refreshGUI() {
        // TODO: if (show) { ?
        treeViewsMap.get(currentPanel).repaint();
    }

    private void createToolBar() {
        toolBar.add(new JToolBar.Separator());
        group = new ButtonGroup();
        for (final TreeView view : treeViews) {
            JToggleButton button = view.getToolBarButton();
            button.setEnabled(false);
            // Move?
            button.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent actionEvent) {
                    CardLayout layout = (CardLayout) treePanel.getLayout();
                    String identifier = view.getIdentifier();

                    currentPanel = identifier;
                    layout.show(treePanel, identifier);
                }
            });
            group.add(button);
            toolBar.add(button);
        }
        // The first (non-separator) is selected by default.
        if (toolBar.size() > 1) {
            ((JToggleButton) toolBar.get(1)).setSelected(true);
        }
    }

    private void addTreeView(TreeView view) {
        String identifier = view.getIdentifier();
        if (treeViewsMap.containsKey(identifier)) {
            throw new RuntimeException("Duplicate exception:" +
                    " You need to change the identifier.");
        }
        if (treeViews.size() == 0) { // Makes the first view the current view.
            currentPanel = identifier;
        }
        treeViewsMap.put(identifier, view);
        treeViews.add(view);
    }


}
