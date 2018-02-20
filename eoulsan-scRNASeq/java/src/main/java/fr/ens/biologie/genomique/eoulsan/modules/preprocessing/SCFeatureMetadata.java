package fr.ens.biologie.genomique.eoulsan.modules.preprocessing;

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
    private int start;
    private int end;
    private String chromosome;
    private String gene = null;
    private String transcript = null;

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
     * Get feature start position
     *
     * @return start position(int)
     */
    public final int getStart() {
        return this.start;
    }

    /**
     * Get feature end position
     *
     * @return end position(int)
     */
    public final int getEnd() {
        return this.end;
    }

    /**
     * Get feature chromosome position
     *
     * @return chromosome (String)
     */
    public final String getChromosome() {
        return this.chromosome;
    }

    /**
     * Get feature transcript ID
     *
     * @return transcript (String)
     */
    public final String getTranscript() {
        return this.transcript;
    }

    /**
     * Test if feature has gene characteristics
     *
     * @return true if feature is a gene, otherwise false
     */
    public final boolean isGene() {
        return (this.getGenomicType().toLowerCase().equals("gene"));
    }

    /**
     * Test if feature is a sparse one
     *
     * @return true if feature is an exon or a CDS
     */
    public final boolean isSparse() {
        return (this.getGenomicType().toLowerCase().equals("exon") || this
            .getGenomicType().toUpperCase().equals("CDS"));
    }

    //
    // Setters
    //

    /**
     * Set genomic type
     *
     * @param type genomic type of the feature
     */
    protected final void setGenomicType(final String type) {
        this.genomicType = type;
    }

    /**
     * Set id of the feature
     *
     * @param id the id of the feature
     */
    protected final void setId(final String id) {
        this.id = id;
    }

    /**
     * Set type of the feature
     *
     * @param type the type of the feature
     */
    protected final void setType(final String type) {
        this.type = type;
    }

    /**
     * Set length of the feature
     *
     * @param length an integer for length of the feature
     */
    protected final void setLength(final int length) {
        this.length = length;
    }

    /**
     * Set gene of the feature
     *
     * @param gene the gene Id of the gene corresponding to the feature
     */
    protected final void setGene(final String gene) {
        this.gene = gene;
    }

    /**
     * Set start position of the feature
     *
     * @param start the start position corresponding to the feature
     */
    protected final void setStart(final int start) {
        this.start = start;
    }

    /**
     * Set end position of the feature
     *
     * @param end the start position corresponding to the feature
     */
    protected final void setEnd(final int end) {
        this.end = end;
    }

    /**
     * Set chromosome position of the feature
     *
     * @param chromosome the chromosome position corresponding to the feature
     */
    protected final void setChromosome(final String chromosome) {
        this.chromosome = chromosome;
    }

    /**
     * Set trabscript of the feature
     *
     * @param transcript the transcript Id of the transcript corresponding to the feature
     */
    protected final void setTranscript(final String transcript) {
        this.transcript = transcript;
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
        this.start = 0;
        this.end = 0;
        this.chromosome = "";
        this.transcript = null;
    }

    /**
     * Write features parameter using a BufferedWriter
     *
     * @param out the Writer
     * @throws IOException error with output file
     */
    public final void printMetadata(BufferedWriter out) throws IOException {

        // Set output for gene
        String output = this.id + "\t" + this.type + "\t" + this.length;
        if (this.isGene()) {
            out.write(output);
            return;
        }

        // Set output for transcripts
        output = output + "\t" + this.gene;

        // Treat sparse features
        if (this.isSparse()) {
            out.write(output + "\t" + this.transcript);
            return;
        }

        out.write(
            this.id + "\t" + this.type + "\t" + this.length + "\t" + this.gene);
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
