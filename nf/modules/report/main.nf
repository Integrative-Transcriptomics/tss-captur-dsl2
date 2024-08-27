/*
    Generates a HTML report using Flask by copying static files and templates needed for the report creation script
    Returns the report in the specified output path
*/
process CREATEREPORT {
    label params.reportLabel
    publishDir "$outputPath/", pattern: "*.tsv", mode: params.pubDirMode

    input:
    val allFigures
    val motifDone
    val outputPath
	
	output:
	path "*.tsv" 

    script:
    """
        cp -r $params.staticfiles static
        cp -r $params.templates templates
        python3 $params.pyCreateReport --path $outputPath --model $params.cnitModel --motifs $params.motifNumber
    """
}

process ZIPREPORT {
    label params.reportLabel
    publishDir "$outputPath/", pattern: "*.zip", mode: params.pubDirMode

    input:
    val report_ready
    val outputPath

	
	output:
	path "*.zip" 

    script:
    """
        timestamp="\$(date +%s)"
        mkdir TSS-CAPTUR-Report-\$timestamp
        cp -r $outputPath/* TSS-CAPTUR-Report-\$timestamp
        zip -r report.zip TSS-CAPTUR-Report-\$timestamp
    """
}