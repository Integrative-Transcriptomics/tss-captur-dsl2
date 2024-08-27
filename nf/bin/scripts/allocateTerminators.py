"""Allocation of terminators for each transcript

This script allocates for each transcript a predicted terminators by 
TransTermHP and RhoTermPredict. 

It recieves the results of both mentioned tools, a CRD file with the 
predicted coordinates for each transcript and a list of all TSS analyzed. 
The goal of the latter file is to avoid the allocating terminators that are
within another transcript. 

This script requires that `pandas`. 

The following functions are found in this script:

    * evaluate_terminator - scores a terminator wrt to the given transcript
    * find_terminators: allocates all terminators that could be allocated
    for each transcript.
    * bin_column: bins the values of the column of RhoTermPredict
    * separate_ids: adapts the identifier for a better parsing
    * get_distances: parses the file containing the maximal extraction length 
    per transcript
    
"""

import re
import pandas as pd
import numpy as np
import argparse


def evaluate_terminator(row, start, end, length_transcript, l_searchspace):
    """Evaluates each terminator wrt. the transcript

    Parameters
    ----------
    row : df_row
        A row of the dataframe that contains the information of the possible
        terminators for this transcript
    start, end: int
        Start and end of the transcript
    length_transcript: int
        The length of transcript
    l_searchspace: int
        The length of the allowed extraction region

    Returns
    -------
    score: float
        The score for the analyzed terminator wrt. to the current transcript
    """


    if row["strand"] == "+":
        transcript_start = start
        transcript_end = end
        terminator_start = row["start"]
    else:
        transcript_start = end
        transcript_end = start
        terminator_start = row["end"]
    # This distance expects the terminator to be outside the predicted area
    distance_to_start = min(
        (abs(terminator_start - transcript_start))/length_transcript, 1)
    # This distance expects the transcript to be close to the close to the end of the transcript
    distance_predicted_end = 1 - \
        (abs(terminator_start - transcript_end)/l_searchspace)
    score = row["binned_score"]/2 + \
        distance_predicted_end/4 + distance_to_start/4

    return score


def find_terminators(row_transcript, df, max_distance):
    """Finds the best terminator for each transcript

    Parameters
    ----------
    row_transcript : df_row
        A row of the dataframe that contains the information of all transcripts
    df: df
        A DataFrame containing all terminators found by TransTermHP and RhoTermPredict
    max_distance: dict
        Dictionary containing the maximal distance used for extraction for each transcript
    
    Returns
    -------
    row_final: df_row
        A modified row_transcript version including the information of the allocated terminator
    """
    row_id, start, end, strand, typeTranscript, utr_end = row_transcript.values
    start, end = int(start), int(end)
    # Depending of the strand, the end or the start of the terminator is taken into account
    is_plus = strand == "+"
    look_for = "start" if is_plus else "end"
    # We still need to compare if the seq we are using is larger than the predicted length provided by QRNA and CNIT, hence the min/max comparison
    start_look_for = min(
        end - max_distance[row_id], start) if not is_plus else max(start, utr_end)
    end_look_for = max(
        start + max_distance[row_id], end) if is_plus else min(end, utr_end)
    # The length of the transcript helps to have the expected distance to the start of the transcript
    length_transcript = abs(start-end)
    length_searched = abs(start_look_for-end_look_for)
    # We recover those found within the same strand and in the interesting search space
    df_subset = df[(df.strand == strand) & (df[look_for].astype("int") >= start_look_for) & (
        df[look_for].astype("int") <= end_look_for)]
    if len(df_subset) == 0:
        return row_transcript
    df_subset["score_term"] = df_subset.apply(
        lambda row: evaluate_terminator(row, start, end, length_transcript, length_searched), axis=1)
    best_terminator_for_transcript = df_subset.sort_values(
        "score_term", ascending=False)[["feature", "attributes", "score_term", "start", "end"]].rename(columns={"attributes": "term_id", "start": "term_start", "end": "term_end"}).head(1).squeeze()
    if best_terminator_for_transcript["feature"] == "RhoTerminator":
        df_hp_term = df_subset[df_subset.feature == "terminator"]
        if is_plus:
            filter_mask_hp = (df_hp_term[look_for] > best_terminator_for_transcript["term_end"]) & (
                df_hp_term[look_for] < best_terminator_for_transcript["term_end"] + 150)
        else:
            filter_mask_hp = (df_hp_term[look_for] < best_terminator_for_transcript["term_start"]) & (
                df_hp_term[look_for] > best_terminator_for_transcript["term_start"] - 150)
        df_hp_term = df_hp_term[filter_mask_hp].sort_values(
            look_for, ascending=is_plus).head(1)
        if len(df_hp_term) != 0:
            best_terminator_for_transcript["feature"] = "RhoTerminatorWithHP"
            best_terminator_for_transcript["additional_term"] = df_hp_term["attributes"].values[0]
            update_coord = "end" if is_plus else "start"
            best_terminator_for_transcript["term_%s" %
                                           update_coord] = df_hp_term[update_coord].values[0]
        else:
            best_terminator_for_transcript["additional_term"] = ""
    else:
        best_terminator_for_transcript["additional_term"] = ""
    row_final = pd.concat(
        [best_terminator_for_transcript, row_transcript])
    return row_final


def bin_column(column):
    """ Bins the column of RhoTermPredict to have the same 
    values as of TransTermHP

    Parameter
    ---------
    column: df_column
        Column containing all confidence values of RhoTermPredict

    Returns
    -------

    binned_column
    
    """
    return pd.cut(
        column, bins=11, labels=np.arange(0, 1.1, 0.1))


def separate_ids(row):
    """Adapts the row identifier

    Parameters
    ----------
    row : df_row
        A row of the dataframe that contains the information of the transcripts


    Returns
    -------
    row: df_row
        A modified row with the adapted identifiers
    """
    string = row[0]
    row["genome"] = re.findall(
        "(\w*_\w*)\S*\|", string)[0]
    row["transcript_id"] = re.findall(
        "\|((?:orphan_|antisense_)\d+)\|", string)[0]
    return row


def get_distances(table):
    """Parses the table containing the maximal extraction length for each transcript 

    Parameters
    ----------
    table: str
        path for the .tsv-file containing the maximal length for each transcript

    Returns
    -------
    df_dict: dict
        dictionary indexed by the genome's ID. Within each key a further dictionary is found. 
        The dictionary contains for each transcript the maximal length 
        that was used for the extraction. 
        
    """
    df_distances = pd.read_csv(table, delimiter="\t", header=None)
    df_distances = df_distances.apply(separate_ids, axis=1)
    df_dict = {k: f.groupby('transcript_id')[1].apply(lambda x: list(
        x)[0]).to_dict() for k, f in df_distances.groupby('genome')}
    return df_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--crd", nargs="+")
    parser.add_argument("--rhoterm", nargs="+")
    parser.add_argument("--nocornac", nargs="+")
    parser.add_argument("--tssAnalyzed")

    args = parser.parse_args()
    args.crd.sort()
    args.rhoterm.sort()
    args.nocornac.sort()
    max_distance = get_distances(args.tssAnalyzed)
    compare_list = zip(args.crd, args.rhoterm, args.nocornac)
    for crd, rhoterm, nocornac in compare_list:
        genome_name = re.split("/", re.split(".crd", crd)[0])[-1]
        rhoterm_df = pd.read_csv(rhoterm, sep="\t", header=None)
        rhoterm_df["binned_score"] = bin_column(rhoterm_df[5])
        nocornac_df = pd.read_csv(nocornac, sep="\t", skiprows=4, header=None)
        nocornac_df["binned_score"] = bin_column(nocornac_df[5])

        terminators_all = pd.concat(
            [rhoterm_df, nocornac_df], ignore_index=True)
        terminators_all.columns = ["seqname", "source", "feature", "start",
                                   "end", "score", "strand", "frame", "attributes", "binned_score"]
        crd_df = pd.read_csv(crd, sep="\t", header=None)
        crd_df.columns = ["transcript_id", "start",
                          "end", "strand", "type", "5UTRend"]
        crd_df = crd_df.apply(
            lambda row: find_terminators(row, terminators_all, max_distance[genome_name]), axis=1, result_type="expand")
        allocated_terminators = crd_df.term_id.unique()
        for term_id in allocated_terminators[pd.notnull(allocated_terminators)]:
            filter_mask = crd_df.term_id == term_id
            filtered_df = crd_df[filter_mask]
            crd_df.loc[filter_mask,
                       "rank_term_score"] = filtered_df.score_term.rank(ascending=False)
        filter_first_rank = crd_df.rank_term_score == 1
        crd_df["start_with_terminator"] = np.where(
            ((crd_df.strand == "-") & filter_first_rank), np.where(crd_df.feature == "RhoTerminator", crd_df["term_start"]-150, crd_df["term_start"]), crd_df["start"])
        crd_df["end_with_terminator"] = np.where(
            ((crd_df.strand == "+") & filter_first_rank),  np.where(crd_df.feature == "RhoTerminator", crd_df["term_end"]+150, crd_df["term_end"]), crd_df["end"])
        # crd_df = crd_df.iloc[:, [6, 9, 7, 8, 0, 1, 2, 5, 10, 3, 4, 11, 12]]
        print(crd_df.columns)
        crd_df = crd_df.loc[:, ["transcript_id", "strand", "type", "start", "end", "feature", "term_id", "score_term", "term_start",
                                "term_end", "rank_term_score", "additional_term", "start_with_terminator", "end_with_terminator"]]
        crd_df["feature"] = crd_df["feature"].replace(
            {"terminator": "Intrinsic", "RhoTerminator": "RhoDep", "RhoTerminatorWithHP": "RhoDepWithHP"})
        crd_df["score_term"] = crd_df["score_term"].round(3)
        crd_df["term_id"] = crd_df["term_id"].str.replace(r'ID=', r'')
        crd_df["additional_term"] = crd_df["additional_term"].str.replace(
            r'ID=', r'')
        crd_df.to_csv("%s_allocated_terminators.tsv" %
                      genome_name, sep="\t", index=False)
