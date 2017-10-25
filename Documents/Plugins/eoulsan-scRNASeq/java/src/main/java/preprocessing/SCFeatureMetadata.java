package preprocessing;

import java.io.BufferedWriter;
import java.io.IOException;

/**
 * This class define an object gathering metadata from annotation file for one feature
 *
 * @author Geoffray Brelurut
 * @since 2017
 */

public class SCFeatureMetadata {

    private String genomicType;
    private String id;
    private String type;
    private int length;
    private String gene = null;

    //
    // Getters
    //

    /**
     * Get feature genomic type
     *
     * @return genomicType
     */
    public final String getGenomicType() {
        return this.genomicType;
    }

    /**
     * Get feature Id
     *
     * @return id
     */
    public final String getId() {
        return this.id;
    }

    /**
     * Get feature type
     *
     * @return type
     */
    public final String getType() {
        return this.type;
    }

    /**
     * Get feature Length
     *
     * @return length
     */
    public final int getLength() {
        return this.length;
    }

    /**
     * Get feature corresponding genId
     *
     * @return gene
     */
    public final String getGene() {
        return this.gene;
    }

    /**
     * Test if feature as gene characteristics
     *
     * @return true if feature is a gene, otherwise false
     */
    public final boolean isGene() {
        return (this.getGenomicType().equals("gene"));
    }

    //
    // Setters
    //

    /**
     * Set genomic type
     *
     * @param type genomic type of the feature
     */
    public final void setGenomicType(final String type) {
        this.genomicType = type;
    }

    /**
     * Set id of the feature
     *
     * @param id the id of the feature
     */
    public final void setId(final String id) {
        this.id = id;
    }

    /**
     * Set type of the feature
     *
     * @param type the type of the feature
     */
    public final void setType(final String type) {
        this.type = type;
    }

    /**
     * Set length of the feature
     *
     * @param length an integer for length of the feature
     */
    public final void setLength(final int length) {
        this.length = length;
    }

    /**
     * Set gene of the feature
     *
     * @param gene the gene Id of the gene corresponding to the feature
     */
    public final void setGene(final String gene) {
        this.gene = gene;
    }

    //
    // Other methods
    //

    /**
     * Clear feature
     */
    public final void clear() {
        this.genomicType = "";
        this.id = "";
        this.type = "";
        this.length = 0;
        this.gene = null;
    }

    /**
     * Write features parameter using a BufferedWriter
     *
     * @param out the Writer
     * @throws IOException error with output file
     */
    public final void printMetadata(BufferedWriter out) throws IOException {
        if (this.isGene()) {
            out.write(this.id + "\t" + this.type + "\t" + this.length);
            return;
        }
        out.write(this.id + "\t" + this.type + "\t" + this.length + "\t" + this.gene);
    }

    //
    // Constructor
    //

    /**
     * Public constructor parameterless
     */
    public SCFeatureMetadata() {
        clear();
    }

    /**
     * Public constructor
     *
     * @param genomicType the genomic type of the Feature
     */
    public SCFeatureMetadata(String genomicType) {
        clear();
        this.genomicType = genomicType;
    }
}
