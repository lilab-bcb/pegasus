#!/usr/bin/env python

from sys import exit

import json
import pkg_resources

from typing import List, Dict, Tuple

import logging
logger = logging.getLogger(__name__)



def load_json_index(input_file: str) -> Dict[str, List[str]]:
    with open(input_file) as fin:
        arr = json.load(fin)
    index_dict = {}
    for value in arr:
        index_dict[value[0]] = value[1]
    return index_dict


def load_chromium_indexes() -> Tuple[dict, dict]:
    # Load chromium index sets
    GA_indexes = load_json_index(pkg_resources.resource_filename("pegasus.check_sample_indexes", "chromium-shared-sample-indexes-plate.json"))
    NA_indexes = load_json_index(pkg_resources.resource_filename("pegasus.check_sample_indexes", "Chromium-i7-Multiplex-Kit-N-Set-A-sample-indexes-plate.json"))
    return GA_indexes, NA_indexes


def load_index_file(index_file: str, GA_indexes: Dict[str, List[str]], NA_indexes: Dict[str, List[str]]) -> List[str]:
    # Load index file
    index_arr = []
    with open(index_file) as fin:
        for line in fin:
            index = line.strip().split(',')[0]
            if index in GA_indexes:
                index_arr.extend([(x, index) for x in GA_indexes[index]])
            elif index in NA_indexes:
                index_arr.extend([(x, index) for x in NA_indexes[index]])
            else:
                index_arr.append((index, 'orig'))
    return index_arr


def hamming(index1: str, index2: str) -> int:
    return sum([index1[i] != index2[i] for i in range(min(len(index1), len(index2)))])


def calc_min_hamming_dist(index_arr: List[str]) -> Tuple[int, int, int]:
    s = len(index_arr)
    min_hd = 1000
    min_i = min_j = -1
    for i in range(s - 1):
        for j in range(i + 1, s):
            hdist = hamming(index_arr[i][0], index_arr[j][0])
            if min_hd > hdist:
                min_hd = hdist
                min_i = i
                min_j = j
    return min_hd, min_i, min_j


def check_index_set(index_set, index_arr, n_mis=1):
    min_hamming = n_mis * 2 + 1

    for sequence in index_set:
        for idx in index_arr:
            hd = hamming(sequence, idx[0])
            if hd < min_hamming:
                return False

    return True


def run_check_sample_indexes(index_file, n_mis=1, n_report=-1):
    GA_indexes, NA_indexes = load_chromium_indexes()
    index_arr = load_index_file(index_file, GA_indexes, NA_indexes)
    min_hd, min_i, min_j = calc_min_hamming_dist(index_arr)

    n_mismatch = (min_hd - 1) // 2
    barcode1 = index_arr[min_i][0] if index_arr[min_i][1] == 'orig' else f"{index_arr[min_i][1]}({index_arr[min_i][0]})"
    barcode2 = index_arr[min_j][0] if index_arr[min_j][1] == 'orig' else f"{index_arr[min_j][1]}({index_arr[min_j][0]})"
    
    logger.info(f"Minimum hamming distance is {min_hd}, achieved between {barcode1} and {barcode2}. A n_mis = {n_mismatch} can be set.")

   
    if n_mismatch < n_mis:
        logger.error(f"Index collision detected in {index_file} with n_mis = {n_mis}!")
    elif n_report > 0:
        logger.info(f"The following 10x RNA indexes are compatible with {index_file}:")
        n_success = 0
        for key, index_set in GA_indexes.items():
            if check_index_set(index_set, index_arr, n_mis):
                print(key)
                n_success += 1
                if n_success >= n_report:
                    break
