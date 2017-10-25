package scDifferentialExpression;

import fr.ens.biologie.genomique.eoulsan.core.InputPorts;
import fr.ens.biologie.genomique.eoulsan.data.DataFormat;
import fr.ens.biologie.genomique.eoulsan.data.DataFormatRegistry;
import fr.ens.biologie.genomique.eoulsan.modules.AbstractModule;



import java.util.HashSet;
import java.util.Set;

/**
 * Created by brelurut on 17/01/17.
 */
public abstract class AbstractSCDifferentialExpression extends AbstractModule {


    /**
     * Parameters Names
     */
    protected static final String CORES = "n.cores";

    /**
     * Default Parameters
     */
    protected static final int DEFAULT_CORES = 1;
    protected int nCores = DEFAULT_CORES;


    //
    // Getters
    //

    /**
     * Get the number of cores to use.
     * @return Returns the number of cores(int)
     */
    protected int getNbCores() { return nCores;}

    //
    // DataFormats
    //

    protected static final DataFormat RESULT_FORMAT =
        DataFormatRegistry.getInstance().getDataFormatFromName("diffexp_result_tsv");

    protected static final DataFormat EXP_FORMAT =
        DataFormatRegistry.getInstance().getDataFormatFromName("filtered_expression_matrix_tsv");

    protected static final DataFormat GENES_FORMAT =
        DataFormatRegistry.getInstance().getDataFormatFromName("genes_metadata_tsv");

    protected DataFormat CELLS_FORMAT =
        DataFormatRegistry.getInstance().getDataFormatFromName("normalized_cells_metadata_tsv");

    //
    // Module Methods
    //

    public abstract InputPorts getInputPorts();
}
