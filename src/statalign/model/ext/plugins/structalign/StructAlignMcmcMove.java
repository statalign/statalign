package statalign.model.ext.plugins.structalign;

import statalign.model.ext.ModelExtension;
import statalign.model.ext.McmcMove;
import statalign.model.ext.plugins.StructAlign;

public abstract class StructAlignMcmcMove extends McmcMove {

	StructAlign owner;
	@Override
	public ModelExtension getOwner() {
		return owner;
	}
}
