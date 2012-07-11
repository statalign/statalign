package statalign.postprocess.plugins.contree.test;

import statalign.postprocess.plugins.TreeNode;
import statalign.postprocess.plugins.contree.CTMain;
import statalign.postprocess.utils.NewickParser;

public class EdgeLengthTester {

    @SuppressWarnings("unused")
	public static void main(String[] args) {
        String[] trees = {
                "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;",
                "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.7)F;"
        };

        NewickParser parser = new NewickParser(trees[0]);
        TreeNode node = parser.parse();

        NewickParser parser2 = new NewickParser(trees[1]);
        TreeNode node2 = parser2.parse();

        CTMain main = new CTMain();
        main.initialize(node, 2);
        main.addNewTree(node);
        main.addNewTree(node2);
        System.out.println(main.constructMajorityTree());
        
    }

}
