package statalign.postprocess.plugins.contree;

import java.util.ArrayList;
import java.util.HashMap;


public class TaxaMap extends HashMap<String, Integer> {

    // Variables

    private ArrayList<String> taxaNames;

    // Functions

    public TaxaMap(int i) {
        super(i);
        taxaNames = new ArrayList<String>(i);
    }

    public TaxaMap() {
        taxaNames = new ArrayList<String>();
    }

    public String getName(int index) {
        return taxaNames.get(index);
    }

    @Override
    public Integer put(String s, Integer integer) {
        Integer i = super.put(s, integer);
        taxaNames.add(s);
        return i;
    }

}
