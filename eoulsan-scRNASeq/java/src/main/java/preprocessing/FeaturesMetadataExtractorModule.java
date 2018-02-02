package preprocessing;

import static fr.ens.biologie.genomique.eoulsan.EoulsanLogger.getLogger;
import static fr.ens.biologie.genomique.eoulsan.core.Modules.renamedParameter;


import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import fr.ens.biologie.genomique.eoulsan.AbstractEoulsanRuntime.EoulsanExecMode;
import fr.ens.biologie.genomique.eoulsan.EoulsanException;
import fr.ens.biologie.genomique.eoulsan.Globals;
import fr.ens.biologie.genomique.eoulsan.annotations.LocalOnly;
import fr.ens.biologie.genomique.eoulsan.bio.io.GFFReader;
import fr.ens.biologie.genomique.eoulsan.bio.io.GTFReader;
import fr.ens.biologie.genomique.eoulsan.bio.GFFEntry;
import fr.ens.biologie.genomique.eoulsan.core.InputPorts;
import fr.ens.biologie.genomique.eoulsan.core.InputPortsBuilder;
import fr.ens.biologie.genomique.eoulsan.core.Modules;
import fr.ens.biologie.genomique.eoulsan.core.OutputPorts;
import fr.ens.biologie.genomique.eoulsan.core.Parameter;
import fr.ens.biologie.genomique.eoulsan.core.StepConfigurationContext;
import fr.ens.biologie.genomique.eoulsan.core.TaskContext;
import fr.ens.biologie.genomique.eoulsan.core.TaskResult;
import fr.ens.biologie.genomique.eoulsan.core.TaskStatus;
import fr.ens.biologie.genomique.eoulsan.core.Version;
import fr.ens.biologie.genomique.eoulsan.data.*;
import fr.ens.biologie.genomique.eoulsan.modules.AbstractModule;
import fr.ens.biologie.genomique.eoulsan.modules.CheckerModule;

import static fr.ens.biologie.genomique.eoulsan.core.OutputPortsBuilder.singleOutputPort;

/**
 * This class define a module that extract information from the annotation file
 *
 * @author Geoffray Brelurut
 * @since 2017
 */

@LocalOnly
public class FeaturesMetadataExtractorModule extends AbstractModule {
    /**
     * Module Name
     */
    private static final String MODULE_NAME = "featuresmetadataextractor";

    /**
     * Parameters Names
     */
    private static final String GENOMIC_TYPE_PARAMETER_NAME = "genomic.type";
    private static final String ATTRIBUTE_ID_PARAMETER_NAME = "attribute.id";

    private static final String OLD_GENOMIC_TYPE_PARAMETER_NAME = "genomictype";
    private static final String OLD_ATTRIBUTE_ID_PARAMETER_NAME = "attributeid";

    private static final String FEATURES_FILE_FORMAT = "features.file.format";

    private static final String MITOCHONDRIAL_TAG_PARAMETER_NAME = "mitochondrial.tag";
    private static final String SPIKE_TAG_PARAMETER_NAME = "spike.tag";

    /**
     * Default Parameters
     */
    private static final String DEFAULT_GENOMIC_TYPE = "exon";
    private static final String DEFAULT_ATTRIBUTE_ID = "Parent";

    private String genomicType = DEFAULT_GENOMIC_TYPE;
    private String attributeId = DEFAULT_ATTRIBUTE_ID;

    private boolean gtfFormat;

    private static final String DEFAULT_MITOCHONDRIAL_TAG = "MT";
    private static final String DEFAULT_SPIKE_TAG = "ERCC";

    private String mitochondrialTag = DEFAULT_MITOCHONDRIAL_TAG;
    private String spikeTag = DEFAULT_SPIKE_TAG;

    // DataFormats
    private static final DataFormat GENES_METADATA_TSV =
            DataFormatRegistry.getInstance().getDataFormatFromName("genes_metadata_tsv");

    //
    // Getters
    //

    /**
     * Get the genomic type.
     *
     * @return Returns the genomicType (String)
     */
    protected String getGenomicType() {
        return this.genomicType;
    }

    /**
     * Get the attribute id.
     *
     * @return Returns the attribute id (String)
     */
    protected String getAttributeId() {
        return this.attributeId;
    }

    /**
     * Test if the annotation file format is GTF
     *
     * @return true if the annotation file format is GTF
     */
    protected boolean isGTFFormat() {
        return this.gtfFormat;
    }

    /**
     * Get the tag used to identify mitochondrial features
     *
     * @return the tag for mitochondrial features
     */
    protected String getMitochondrialTag() {
        return this.mitochondrialTag;
    }

    /**
     * Get the tag used to identify spike-ins' features
     *
     * @return the tag for spike-ins' features
     */
    protected String getSpikeTag() {
        return this.spikeTag;
    }

    //
    // Module Methods
    //

    @Override
    public String getName() {
        return MODULE_NAME;
    }

    @Override
    public Version getVersion() {
        return Globals.APP_VERSION;
    }

    @Override
    public InputPorts getInputPorts() {

        final InputPortsBuilder builder = new InputPortsBuilder();

        builder.addPort("featuresannotation", this.gtfFormat ?
                DataFormats.ANNOTATION_GTF : DataFormats.ANNOTATION_GFF);

        return builder.create();
    }

    @Override
    public OutputPorts getOutputPorts() {
        return singleOutputPort("genesoutput", GENES_METADATA_TSV);
    }

    @Override
    public void configure(final StepConfigurationContext context,
                          final Set<Parameter> stepParameters) throws EoulsanException {

        for (Parameter p : stepParameters) {

            // Check if parameter is deprecated

            switch (p.getName()) {

                case OLD_GENOMIC_TYPE_PARAMETER_NAME:
                    renamedParameter(context, p, GENOMIC_TYPE_PARAMETER_NAME);
                case GENOMIC_TYPE_PARAMETER_NAME:
                    this.genomicType = p.getStringValue();
                    break;

                case OLD_ATTRIBUTE_ID_PARAMETER_NAME:
                    renamedParameter(context, p, ATTRIBUTE_ID_PARAMETER_NAME);
                case ATTRIBUTE_ID_PARAMETER_NAME:
                    this.attributeId = p.getStringValue();
                    break;

                case MITOCHONDRIAL_TAG_PARAMETER_NAME:
                    this.mitochondrialTag = p.getStringValue();
                    break;

                case SPIKE_TAG_PARAMETER_NAME:
                    this.spikeTag = p.getStringValue();
                    break;

                case FEATURES_FILE_FORMAT:
                    switch (p.getLowerStringValue()) {

                        case "gtf":
                            this.gtfFormat = true;
                            break;

                        case "gff":
                            this.gtfFormat = false;
                            break;
                        case "gff3":
                            this.gtfFormat = false;
                            break;

                        default:
                            Modules.badParameterValue(context, p, "Unknown annotation file format");
                            break;
                    }
                    break;

                default:
                    Modules.unknownParameter(context, p);
            }
        }

        // Log Step parameters
        getLogger().info("In " + getName() + ", mitochondrial tag=" +
                this.mitochondrialTag + ", spike tag=" + this.spikeTag);
        getLogger().info("In " + getName() + ", genomic type=" +
                this.genomicType + ", attribute ID=" + this.attributeId);
    }


    @Override
    public TaskResult execute(final TaskContext context,
                              final TaskStatus status) {
        try {
            final Data featuresAnnotationData =
                    context.getInputData(isGTFFormat() ? DataFormats.ANNOTATION_GTF : DataFormats.ANNOTATION_GFF);
            final Data geneMetadata = context.getOutputData(GENES_METADATA_TSV, featuresAnnotationData);

            // Get annotation file
            final DataFile annotationFile = featuresAnnotationData.getDataFile();

            // Get final metadata file
            final DataFile metadataFile = geneMetadata.getDataFile();

            // Write metadata file
            extractMetadata(annotationFile, this.mitochondrialTag, this.spikeTag,
                    this.genomicType, this.attributeId, this.gtfFormat, metadataFile);

            // Write log file
            return status.createTaskResult();
        } catch (EoulsanException e) {
            return status.createTaskResult(e, e.getMessage());
        } catch (FileNotFoundException e) {
            return status.createTaskResult(e, "File not found: " + e.getMessage());
        } catch (IOException e) {
            return status.createTaskResult(e, "Error while reading annotation file: " + e.getMessage());
        }
    }

    /**
     * Extract and Write features metadata from annotation file
     *
     * @param annotations the annotation file
     * @param mtTag       tag for mitochondrial features
     * @param spikeTag    tag for spike in features
     * @param genomicType genomic type of the features to consider
     * @param attributeId Id to design the feature
     * @param gtfFormat   boolean indicating if the annotation file is gtf or not
     * @param outFile     output file to write data
     * @throws IOException      if encouters problem with input or output file
     * @throws EoulsanException if features as no Id
     */
    private void extractMetadata(DataFile annotations, final String mtTag, final String spikeTag,
                                 final String genomicType, final String attributeId, final boolean gtfFormat,
                                 DataFile outFile) throws IOException, EoulsanException {

        try (final GFFReader annotationReader =
                     gtfFormat ? new GTFReader(annotations.open()) : new GFFReader(annotations.open());
             BufferedWriter out = new BufferedWriter(new FileWriter(outFile.toFile()))) {

            // Initialise gathering variables
            Map<String, SCFeatureMetadata> features = new HashMap<>();
            Map<String, String> parents = new HashMap<>();

            //Read annotation file
            for (final GFFEntry anno : annotationReader) {

                // Treat non genic case, saving gene ID
                if (!genomicType.equals("gene")) {
                    parents = saveGeneId(anno, parents, gtfFormat);
                }

                // Extract information
                if (genomicType.equals(anno.getType())) {

                    String type;

                    // get Id
                    final String featureId = anno.getAttributeValue(attributeId);

                    if (featureId == null) {
                        throw new EoulsanException("Feature " + genomicType + " does not contain a "
                                + attributeId + " attribute");
                    }

                    // get length
                    final int start = anno.getStart();
                    final int end = anno.getEnd();
                    final int length = end - start + 1;

                    // get feature type
                    if (anno.getSeqId().equals(mtTag)) {
                        type = "mitochondrial";
                    } else if (anno.getSeqId().equals(spikeTag)) {
                        type = "spike";
                    } else {
                        type = "other";
                    }

                    // Test if object is referenced (sparse features) else create it
                    if (!features.containsKey(featureId)) {
                        SCFeatureMetadata feature = new SCFeatureMetadata(genomicType);
                        feature.setLength(length);
                        feature.setType(type);
                        feature.setId(featureId);

                        features.put(featureId, feature);
                    } else {
                        SCFeatureMetadata feature = features.get(featureId);
                        final int oldLength = feature.getLength();
                        feature.setLength(oldLength + length);
                        features.put(featureId, feature);
                    }
                }
            }
            for (String id : features.keySet()) {
                SCFeatureMetadata feature = features.get(id);
                if (!feature.isGene()) {
                    feature.setGene(parents.get(id));
                }
            }

            // Write headers in output file
            out.write( "Type" + "\t" + "Length");
            if (!genomicType.equals("gene")) {
                out.write("\t" + "GeneID");
            }
            // Write metadata in output file
            out.newLine();
            for (String id : features.keySet()) {
                features.get(id).printMetadata(out);
                out.newLine();
            }
        }
    }

    /**
     * Save genes id in a already existing Map
     *
     * @param anno      a GFFentry object
     * @param map       Map to save data
     * @param gtfFormat boolean indicating whether annotations are from a gtf file or not
     * @return Map completed with key : feature id, value : corresponding gene id
     */
    private Map<String, String> saveGeneId(GFFEntry anno, Map<String, String> map, boolean gtfFormat) {
        if (gtfFormat && anno.getType().matches("(.)*transcript(.)*")) {
            map.put(anno.getAttributeValue("transcript_id"), anno.getAttributeValue("gene_id"));
            return map;
        }
        if (anno.getType().matches("(.)*transcript(.)*")) {
            map.put(anno.getAttributeValue("ID"), anno.getAttributeValue("Parent"));
            return map;
        }
        return map;
    }
}
