from flask import Flask, render_template
import argparse
import os
import re
import pandas as pd
import numpy as np
from flask_frozen import Freezer
from functools import reduce

app = Flask(__name__)
freezer = Freezer(app)


toFillNA = "NA"


def createRows(row):
    row_id_all = row[0]
    row["transcript_id"] = re.findall(
        "\|((?:orphan_|antisense_)\d+)\|", row_id_all)[0]
    row["position"] = int(re.findall("\|Start:(\d+)\|", row_id_all)[0])
    row["strand"] = re.findall("\|Strand:(\+|\-)", row_id_all)[0]
    row["transcript_length"] = row[1]
    return row.drop([0, 1])


def modifyMFE(row):
    row["transcript_id"] = row[0].replace(">", "")
    row["mfe"] = row[3]
    return row[["transcript_id", "mfe"]]


def createShort(overview, ignored):
    ignoredCount = len(ignored)
    dict_return = reduce(lambda d, x: update_dict(d, x[0], x[4]), overview, {
                         "antisense": {"COD": 0, "RNA": 0}, "orphan": {"COD": 0, "RNA": 0}})
    dict_return.update({"Ignored": ignoredCount})

    return dict_return


def update_dict(d, type_tss, class_tss):
    temp_type = type_tss.split("_")[0]
    current_subdict = d.get(temp_type, {})
    current_count = current_subdict.get(class_tss, 0)
    current_subdict.update({class_tss: current_count+1})
    d.update({temp_type: current_subdict})
    return d


def createOverviewData(pathToData, genome):
    dir_tss_analyzed = os.path.join(pathToData, "Queries", "tss_analyzed.tsv")
    dir_final_class = os.path.join(
        pathToData, "Classification", "%s_final_classification.tsv" % genome)
    dir_terminators = os.path.join(
        pathToData, "Terminators", "%s_allocated_terminators.tsv" % genome)
    dir_sec_struct = os.path.join(
        pathToData, "SecondaryStructure", genome, "%s.tsv" % genome)
    dir_motif_analysis = os.path.join(pathToData, "MotifAnalysis", genome)
    dir_motif_analysis_seqs = os.path.join(
        dir_motif_analysis, "%s_seqs_meme_results.tsv" % genome)
    dir_motif_summaries = os.path.join(
        dir_motif_analysis, "%s_motifs_summary.tsv" % genome)
    tss_analyzed = pd.read_csv(dir_tss_analyzed, sep="\t", header=None)
    tss_analyzed = tss_analyzed[tss_analyzed[0].str.contains(
        genome)].apply(createRows, axis=1).set_index(["transcript_id"])
    classification = pd.read_csv(
        dir_final_class, sep="\t", index_col="transcript_id")[["winner", "PredictionStart", "PredictionEnd"]]
    terminators = pd.read_csv(
        dir_terminators, sep="\t", index_col="transcript_id")[["feature", "term_start", "term_end", "start_with_terminator", "end_with_terminator"]]
    mfe = pd.read_csv(dir_sec_struct, sep="\t",
                      header=None).apply(modifyMFE, axis=1).set_index("transcript_id")
    motifs = pd.read_csv(dir_motif_analysis_seqs,
                         sep="\t", index_col="transcript_id")[["Motif ID", "Motif Start", "p-value"]]
    summary_motifs = pd.read_csv(dir_motif_summaries,
                                 sep="\t")[["Motif ID", "Motif Name", "Evalue", "Occurrences"]]
    motifs = motifs.groupby(motifs.index).agg(list)
    motifs["summary"] = motifs.apply(lambda x: [list(a) for a in zip(
        x["Motif ID"], x["Motif Start"], x["p-value"])], axis=1)
    motifs.drop(["Motif ID", "Motif Start", "p-value"], axis=1, inplace=True)
    all_dfs = pd.concat(
        [tss_analyzed, classification, terminators, mfe, motifs], join="outer", axis=1)
    all_dfs["GeneStart"] = all_dfs.start_with_terminator.where(
        (all_dfs.strand == "-") | ((all_dfs.strand == "+") & (all_dfs.PredictionStart > all_dfs.end_with_terminator)), all_dfs["PredictionStart"])
    all_dfs["GeneEnd"] = all_dfs.end_with_terminator.where(
        (all_dfs.strand == "+") | ((all_dfs.strand == "-") & (all_dfs.PredictionEnd < all_dfs.start_with_terminator)), all_dfs["PredictionEnd"])
    all_dfs["GeneLength"] = abs(all_dfs.GeneEnd - all_dfs.GeneStart)+1
    temp = all_dfs[["GeneEnd", "position", "GeneStart"]]
    temp_max = temp.max(axis=1)
    temp_min = temp.min(axis=1)
    all_dfs["GeneLengthWithUTR"] = temp_max - temp_min+1
    all_dfs["mfe"] = all_dfs.mfe.astype("float")
    all_dfs["url_for_image"] = np.where(
        (all_dfs.winner == "RNA") | (all_dfs.winner == "COD"), all_dfs.index.astype("str") + "_ss.jpg", "")
    all_dfs = all_dfs.drop(
        ["start_with_terminator", "end_with_terminator", "PredictionEnd", "PredictionStart"], axis=1).reset_index().fillna(toFillNA)
    all_dfs.to_csv("%s_overview.tsv" % genome, sep="\t", index=False)
    return [all_dfs.values.tolist(), summary_motifs.values.tolist()]


@ app.context_processor
def give_genomes():
    overviewData = data["overview"]
    classifiedData = data["classification"]
    avoided = data["avoided"]
    terminators = data["terminators"]
    summaryMotifs = data["summary_motifs"]
    short_description = data["short_description"]
    return dict(genomes=genomes, outputPath=outputPath, short_description=short_description, avoidedTSS=avoided, classified=classifiedData, terminators=terminators, overviewData=overviewData, summaryMotifs=summaryMotifs)


@ app.route('/overview.html')
def overview():
    return render_template('overview.html')


def getClassification(outputPath, genome):
    dir_final_class = os.path.join(
        outputPath, "Classification", "%s_final_classification.tsv" % genome)
    data = pd.read_csv(
        dir_final_class, sep="\t").fillna(toFillNA).values.tolist()
    return data


@ app.route('/classified.html')
def classified():
    return render_template('classified.html')


def getTerminators(outputPath, genome):
    dir_terminators = os.path.join(
        outputPath, "Terminators", "%s_allocated_terminators.tsv" % genome)
    terminators = pd.read_csv(
        dir_terminators, sep="\t").fillna(toFillNA).values.tolist()
    return terminators


@ app.route('/terminators.html')
def terminators():
    return render_template('terminators.html')


def getAvoidedTSS(outputPath, genome):
    dir_avoided = os.path.join(
        outputPath, "Queries", "tss_avoided_in_%s.tsv" % genome)
    try:
        df_avoided = pd.read_csv(dir_avoided, sep="\t", header=None)
        avoided = df_avoided.values.tolist()

    except:
        avoided = pd.DataFrame([]).values.tolist()
    return avoided


@ app.route('/avoided.html')
def avoided():
    return render_template('avoided.html')


@ app.after_request
def add_header(r):
    """
    Add headers to both force latest IE rendering engine or Chrome Frame,
    and also to cache the rendered page for 10 minutes.
    """
    r.headers["Cache-Control"] = "no-cache, no-store, must-revalidate"
    r.headers["Pragma"] = "no-cache"
    r.headers["Expires"] = "0"
    r.headers['Cache-Control'] = 'public, max-age=0'
    return r


if __name__ == '__main__':
    index = 0
    parser = argparse.ArgumentParser()
    parser.add_argument("--path")
    args = parser.parse_args()
    outputPath = args.path
    print("Parsing files in %s" % outputPath)
    genomes_list = os.path.join(args.path, "genomes_text.txt")
    with open(genomes_list) as f:
        genomes = [i.strip() for i in f.readlines()]
    genomes.sort()
    print("Genomes found: %s" % ", ".join(genomes))
    data = {"avoided": {}, "overview": {},
            "classification": {}, "terminators": {}, "summary_motifs": {}, "short_description": {}}
    for it, g in enumerate(genomes):
        print("Parsing data: {} % ".format(100 * float(it) / len(genomes)))
        data["avoided"][g] = getAvoidedTSS(outputPath, g)
        data["terminators"][g] = getTerminators(outputPath, g)
        data["classification"][g] = getClassification(outputPath, g)
        temp_overview = createOverviewData(outputPath, g)
        data["overview"][g] = temp_overview[0]
        data["summary_motifs"][g] = temp_overview[1]
        data["short_description"][g] = createShort(
            temp_overview[0], data["avoided"][g])
    print("Parsing completed. Freezing Interface.")

    genome = genomes[index]
    app.jinja_env.auto_reload = True
    app.config['TEMPLATES_AUTO_RELOAD'] = True
    app.config["FREEZER_DESTINATION"] = os.path.join(
        args.path, "Interface")
    freezer.freeze()
