"""

Author           : Paul Quinn
Email            : paul.quinn-2@postgrad.manchester.ac.uk
Github           : https://github.com/peq123/
Supervisor       :  May Tassabehji, Peter Causey-Freeman
Supervisor Email : m.tassabehji@manchester.ac.uk, peter.causey-freeman@manchester.ac.uk

Description:
This module is used to accumulate human genomic data from public API servies, particularly focusing on variants. 

Usage:



"""
from __future__ import absolute_import
from __future__ import print_function

import os
import sys
import time
import requests  # for web calls
import urllib  # for url encoding
import html  # for escaping > <
import re  # regex
import logging
import json
import pandas as pd
import math

from datetime import datetime

import urllib  # for url encoding
if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO

logger = logging.getLogger(__name__)

logger.setLevel(logging.WARNING)
# logger.setFormatter(logging.Formatter("%(name)s:%(levelname)s: %(message)s"))

# class for genes


class parameters:

    active_assembly = "GRCh38"
    # this allows the user to store the variants locally for fast loading
    store_local = True


# base api class
class api:

    __assemblies = None
    __name = ""
    request_history = None
    history_limit = 10

    def __init__(self, *args):
        raise Exception("This class can not have instances")

    @classmethod
    def access(cls, url, isjson=True, timeout=1000, headers=None, data=None):
        """ call url and retrieve json """
        if cls.request_history is None:
            cls.request_history = []
        try:
            if headers is None:
                headers = {}
            if data is None:
                r = requests.get(url, timeout=timeout, headers=headers)
            else:
                r = requests.post(url,
                                  timeout=timeout,
                                  headers=headers,
                                  data=data)
            cls.request_history = [
                r
            ] + cls.request_history[0:cls.history_limit - 1]
        except Exception as ex:
            cls.request_history = [
                ex
            ] + cls.request_history[0:cls.history_limit - 1]

            raise ex

        if not r.ok:
            if r.status_code == 500:
                raise Exception(
                    "HTTPError 500: Website scripting error, this maybe caused by requesting url. Please check any search terms."
                )
            else:
                r.raise_for_status()

        if isjson:

            return r.json()
        else:

            return r.text

    @classmethod
    def rename_keys(cls, thedict, prefix="", suffix=""):
        if type(thedict) is dict:
            return {
                "{0}{1}{2}".format(prefix, k, suffix): v
                for k, v in thedict.items()
            }
        elif type(thedict) is str:
            return "{0}{1}{2}".format(prefix, thedict, suffix)

    @classmethod
    def genomic_hgvs(cls, cdna_hgvs, genomic_coordinates=None, returntype=str):
        """ Returns hex encoded string to be used as dataframe index. Will translate to used reference assembly using liftover """

        if type(genomic_coordinates) is not dict:
            raise Exception(
                "Function expects position as a dictionary, with at least one entry that has the assembly as the key"
            )

        if type(cdna_hgvs) is not str or cdna_hgvs == "":
            raise Exception(
                "Function expects the cDNA HGVS description for the variant")

        cdna_hgvs = html.unescape(cdna_hgvs)
        changetype = re.search("(?<=[^\:]\:)[a-zA-Z]", cdna_hgvs)

        if changetype is None or changetype[0].lower() not in ["c", "g"]:
            return cdna_hgvs  # fail
        needs_complement = changetype[0].lower() in ["c"]
        completechange = re.search("(?<=[\d|\)])[a-z|A|C|G|T|>]{3,}",
                                   cdna_hgvs)

        if needs_complement:
            # complement
            completechange = completechange[0].translate(
                str.maketrans({
                    'A': 'T',
                    'C': 'G',
                    'G': 'C',
                    'T': 'A'
                }))
        else:
            completechange = completechange[0]

        return

    @classmethod
    def __default_assemblies(cls):
        return pd.DataFrame(
            index=["GRCh38", "GRCh37"],
            data={"names": [["grch38", "hg38"], ["grch37", "hg19", "hg37"]]},
        )

    @classmethod
    def get_assemblies(cls):

        # the __assemblies hasn't been initialised, do so now
        if type(cls.__assemblies) is not pd.DataFrame:
            cls.__assemblies = cls.__default_assemblies()

        return cls.__assemblies

    @classmethod
    def get_active_assembly(cls):
        df = cls.get_assemblies()
        """ Retrieves the active assembly """
        return df[df["names"].apply(lambda x: parameters.active_assembly.lower(
        ) in [i.lower() for i in x])]

    @classmethod
    def add_assembly(cls, index, names, is_active=False):
        """ Adds assembly to the assembly dataframe. Index has to be unique. """

        # the __assemblies hasn't been initialised, do so now
        if type(cls.__assemblies) is not pd.DataFrame:
            cls.__assemblies = cls.__default_assemblies()

        # type checks
        if type(index) is not str:
            raise Exception('"index" is expected to be a string.')
        if type(names) is not list:
            raise Exception('"names" is expected to be a list of values.')
        if type(is_active) is not bool:
            raise Exception('"is_active" is expected to be a boolean.')

        # placeholder for future use
        if is_active:
            logging.warn(
                "Setting current assembly is currently not supported and is a placeholder for future functionality."
            )

        if cls.assembly_index(index)[0] is None:
            raise Exception(
                "Index already exists in the assemblies dataframe.")

        # add new assembly
        cls.__assemblies = pd.concat([
            cls.__assemblies,
            pd.DataFrame(
                index=[index],
                data=[{
                    "names": names,
                    "is_active": False,
                    "is_default": False
                }],
            ),
        ])

    @classmethod
    def remove_assembly(cls, index):
        """ Attempts to remove the specified assembly index from the assembly dataframe. Cannot remove from defaults."""
        if type(cls.__assemblies) is not pd.DataFrame:
            cls.__assemblies = cls.__default_assemblies()

        if type(index) is not str:
            raise Exception('"index" is expected to be a string.')

        # if the index exists this will return the case sensitive index otherwise index will be set to None
        index, is_default = cls.assembly_index(index)

        if is_default:
            raise Exception(
                "This index is part of the default assemblies. Cannot be removed."
            )

        # does it exist?
        if index is None:
            raise Exception(
                "This index does not exist in the assemblies dataframe.")

        # dataframe indices are case sensitive therefore extract the correct case from the dictionary prepared above
        cls.__assemblies.drop(index)

    @classmethod
    def assembly_index(cls, index, default=None):
        """ returns the actual key for the provided index if the assembly exists """
        if type(cls.__assemblies) is not pd.DataFrame:
            cls.__assemblies = cls.__default_assemblies()
        # first attempt to find the index
        possibles = cls.__assemblies[cls.__assemblies["names"].apply(
            lambda x: index.lower() in [i.lower() for i in x])]

        # dataframe indices are case sensitive, create dictionary of indices with the lowercase as keys
        if len(possibles.index) > 0:
            return (possibles.index[0].lower(), parameters.active_assembly in
                    cls.__assemblies.loc[possibles.index[0], "names"])
        else:
            return default, False

    @classmethod
    def translate_assembly_name(cls, obj):

        if type(obj) is str:
            return cls.assembly_index(obj)[0]
        elif type(obj) is list:

            return [cls.assembly_index(i)[0] for i in obj if type(i) is str]

        elif type(obj) is dict:

            return {cls.assembly_index(k)[0]: v for k, v in obj.items()}

    @classmethod
    def __aggregate_frame(cls, frame):
        """ Takes a pandas dataframe and returns an aggregate form. Best used inconjunction with the groupby function. """
        tobereturned = {}

        for c in frame.columns:
            values = []

            thistype = str
            for v in frame[c].values:
                # convert the valuest
                if type(v) is list:
                    values.extend(v)
                else:
                    values.append(str(v))
            theset = set(values)

            if len(theset) > 1:
                # the set has more than one value therefore aggregate
                tobereturned[c] = values
            else:
                tobereturned[c] = values[0]

        return pd.DataFrame([tobereturned])

    @classmethod
    def combine_duplicates(cls, frame, columns):
        """ Uses the pandas dataframe.groupby method to remove duplicates based upon grouping on the specified columns. """
        if type(columns) is not list and type(columns) is not str:
            raise Exception(
                "Columns is expected to be either a string or a list")

        return (frame.groupby(columns).apply(
            cls.__aggregate_frame).reset_index(drop=True))

    @classmethod
    def extract_hgvsc(cls, df, fromname, prefix=""):
        return df
        """ Updates the dataframe to have a "change" column and an "change_issue" column extracted from the provided column  """
        # this removes the reference transcript from the beginning and any ammino acid labels

        if type(df) is not pd.DataFrame:
            raise Exception("Function expects the 'df' as a pandas DataFrame.")

        if type(fromname) is not str:
            raise Exception(
                "Function expects the 'fromname' column name as a string.")

        if type(prefix) is not str:
            prefix = ""

        if not prefix == "" and not prefix[-1] == "_":
            prefix = prefix + "_"  # add underscore to the end

        df["{0}change".format(prefix)] = df.apply(
            lambda x: html.unescape(re.search("c.[^\s]*", x[fromname])[0])
            if ":c." in x[fromname] else "[" + x[fromname] + "]",
            axis=1)
        df["{0}valid_variant".format(prefix)] = df.apply(lambda x: re.search(
            "[\[|\]|?|=|,]", x["{0}change".format(prefix)]) is None,
                                                         axis=1)
        return df

    @classmethod
    def save_local(cls, df, datatype=""):
        """
        This attempts to load the variants from a local source at the script location under the following directory tree:
        DIR / RawData / APIname / type
        There are no checks on this if the user edits these files.
        """
        if type(df) is not pd.DataFrame:
            raise Exception("Function expects a Pandas DataFrame.")

        if type(datatype) is not str:
            raise Exception("Function expects datatype as a string.")

        directory = "{0}{1}RawData{1}{2}".format(
            os.path.dirname(os.path.realpath(__file__)), os.path.sep,
            cls.__name__)
        if not os.path.isdir(directory):
            os.makedirs(directory)

        df.to_csv("{0}{1}{2}.csv".format(directory, os.path.sep, datatype))

    @classmethod
    def load_local(cls, datatype=""):
        logger.debug("{0}{1}RawData{1}{2}{1}{3}.csv".format(
            os.path.dirname(os.path.realpath(__file__)), os.path.sep,
            cls.__name__, datatype))
        logger.debug(
            os.path.exists("{0}{1}RawData{1}{2}{1}{3}.csv".format(
                os.path.dirname(os.path.realpath(__file__)), os.path.sep,
                cls.__name__, datatype)))
        return cls.load_path("{0}{1}RawData{1}{2}{1}{3}.csv".format(
            os.path.dirname(os.path.realpath(__file__)), os.path.sep,
            cls.__name__, datatype))

    @classmethod
    def load_path(cls, path):

        df = pd.read_csv(
            path, index_col=0)  # load the csv file with index from column 0
        # load the data
        return cls.parse_frame(df)  # df.apply(cls.parse_data)

    @classmethod
    def parse_frame(cls, df):

        for c in df.columns:
            if (df[c].apply(lambda x: str(x).startswith("[") or str(x).
                            startswith("{")).any()):
                df[c] = df[c].apply(cls.parse_data)

        return df

    @classmethod
    def parse_data(cls, x):

        strver = str(x)

        # extract out [] objects

        if strver.startswith('"[') or strver.startswith('"{'):
            strver = strver[1:]

        if strver.endswith(']"') or strver.endswith('}"'):
            strver = strver[0:-1:]
        parsedall = False

        if strver.startswith("[") or strver.startswith("{"):
            strver = strver.replace('"{"', '{"')
            strver = strver.replace('"}"', '"}')
            strver = strver.replace("'{'", "{'")
            strver = strver.replace("'}'", "'}")

            # now there could be nested quotes
            strver = strver.replace("\"{'", "{'")
            strver = strver.replace("'}\"", "'}")

            strver = strver.replace(": None", ": ''")
            parts = strver.split('"')
            for i, p in enumerate(parts):
                if parts[i].count("'") == 1:
                    parts[i] = p.replace("'", urllib.parse.quote("'"))
                else:
                    parts[i] = p.replace("'", '"')

            strver = '"'.join(parts)
            try:
                return json.loads(strver)
            except Exception as ex:
                #logger.warn("Issue loading json {0} {1}\t{2}\t{3}".format(len(parts), type(x), x, strver))
                # just return x
                return x
                #raise Exception(ex)
        else:
            return x


class GeneNames(api):
    """ Use this API class to access information provided by the HGNC API """
    @classmethod
    def info(cls):
        if "terms" not in cls.__dict__ or type(cls.terms) is not dict:
            try:
                results = cls.access(
                    "http://rest.genenames.org/info".format(),
                    headers={"Accept": "application/json"},
                )
                cls.terms = results
            except:
                raise Exception(
                    "Unexpected problem occured when trying to retrieve HGNC API information."
                )

    # specific functions
    @classmethod
    def fetch(cls, value, term="symbol"):

        cls.info()

        if type(term) is not str or term == "":
            print(term, type(term))
            raise Exception("Term needs to be a valid string")

        if term not in cls.terms["storedFields"]:
            raise Exception(
                "The term {0} is not a valid stored field to fetch from HGNC.".
                format(term))

        results = cls.access(
            "http://rest.genenames.org/fetch/{0}/{1}".format(term, value),
            headers={"Accept": "application/json"},
        )
        try:
            results = results["response"]
        except:
            raise Exception("Unexpected format from HGNC GeneNames.")

        if "numFound" in results and results["numFound"] == 0:
            raise Exception("Fetching the {0} {0} returned no results.".format(
                term, value))

        return results

    @classmethod
    def search(cls, value, term=""):
        cls.info()

        if term not in cls.terms["searchableFields"] and not term == "":
            raise Exception(
                "The term {0} is not a valid stored field to fetch from HGNC.".
                format(term))

        results = cls.access(
            "http://rest.genenames.org/search/{0}{1}{2}".format(
                term, "/" if not term == "" else "", value),
            headers={"Accept": "application/json"},
        )
        try:
            results = results["response"]
        except:
            raise Exception("Unexpected format from HGNC GeneNames.")

        return results


class VariantValidator(api):
    """ Use this API class to validate the variant """
    @classmethod
    def validator(cls, assembly, description, transcript="all"):
        """ Use API to validate the transcript """

        if type(assembly) is not str or assembly == "":
            raise Exception("Function expects description as string.")

        if type(description) is not str or description == "":
            raise Exception("Function expects description as string.")

        if type(transcript) is list:
            transcript = "|".join(transcript)

        if type(transcript) is not str or transcript == "":
            raise Exception("Function expects description as string.")

        return cls.access(
            "https://rest.variantvalidator.org/VariantValidator/variantvalidator/{0}/{1}/{2}?content-type=application%2Fjson"
            .format(assembly, description, transcript),
            headers={"accept": "application/json"})

    @classmethod
    def gene2transcripts(cls, query):
        """ Provide HGNC gene symbol or transcript to return the refseq transcripts. """

        if type(query) is not str or query == "":
            raise Exception(
                "Function expects the query as either a HGNC gene symbol or transcript as string."
            )

        return cls.access(
            "https://rest.variantvalidator.org/VariantValidator/tools/gene2transcripts/{0}?content-type=application%2Fjson"
            .format(query),
            headers={"accept": "application/json"})

    @classmethod
    def genomic_coordinates(cls, assembly_name, description, return_assembly):

        results = cls.validator(assembly_name, description)

    @classmethod
    def get_transcript_versions(cls, symbol):

        results = cls.gene2transcripts(symbol)
        if "error" in results:
            raise Exception("Unexpected result from API, please check symbol")

        try:
            return {
                i["reference"].split(".")[0]: i["reference"].split(".")[1]
                if "." in i["reference"] else i["reference"]
                for i in results["transcripts"]
            }
        except:
            raise Exception("Unexpected result from API, please check symbol")


class Ensembl(api):

    include_assembly = False
    default_assembly = ""
    versions = None
    apiurl = ""

    @classmethod
    def __apiurl(cls):
        # tests the active assembly to determine the url
        if cls.default_assembly == "":
            try:
                results = cls.access(
                    "https://rest.ensembl.org/info/assembly/homo_sapiens?content-type=application/json"
                )
                cls.default_assembly = results["default_coord_system_version"]
                cls.versions = [
                    i.lower() for i in results["coord_system_versions"]
                ]
            except:
                logger.warn(
                    "Ensembl: Error retrieving current assembly version")
                raise Exception("Error retrieving current assembly version")

        if cls.default_assembly == parameters.active_assembly:
            return "https://rest.ensembl.org/"
        elif parameters.active_assembly.lower() in cls.versions:
            return "https://{0}.rest.ensembl.org/".format(
                parameters.active_assembly.lower())

    @classmethod
    def get_xrefs(cls, id):
        pass
        results = cls.access(
            "{0}/xrefs/id/{1}?content-type=application/json".format(
                cls.__apiurl(), id))
        returninfo = {"symbol": "", "entrez": ""}

        # process the information to extract the needed information
        for ref in results:
            if ref["dbname"] == "HGNC":
                returninfo["symbol"] = ref["display_id"]
                pass
            elif ref["dbname"] == "EntrezGene":
                returninfo["entrez"] = ref["primary_id"]

        return returninfo

    @classmethod
    def get_transcripts(cls, thisgene):
        """ Retrieves list of transcripts + exons from Ensembl """
        if type(thisgene) is not gene:
            raise Exception("Function expects gene object")

        if parameters.store_local:
            try:
                # try and load local copy

                return_var = cls.load_local(
                    str(thisgene.ensembl_id) + "_" +
                    parameters.active_assembly + "_transcripts")
                logger.info("Entrez: Transcript dataset loaded from file.")
                return return_var
            except Exception as ex:
                # if that fails retrieve from api
                pass

        logger.info("Ensembl: Requesting overlap for transcripts and exons")
        results = cls.access(
            "{0}overlap/id/{1}?feature=transcript;feature=exon;content-type=application/json"
            .format(cls.__apiurl(), thisgene.ensembl_id))

        # transcripts = {
        #    i["id"]: transcript(i, "ensembl") for i in results if i["Parent"] == id
        # }

        logger.info("Ensembl: Parsing data")
        datalist = [
            cls.rename_keys(i, prefix="ensembl_") for i in results
            if i["Parent"] == thisgene.ensembl_id
        ]

        for t in datalist:
            t["ensembl_exons"] = {
                i["id"]: i
                for i in results if i["Parent"] == t["ensembl_id"]
            }

        # index is set to the genomic coordinates
        transcripts = pd.DataFrame(
            data=datalist,
            #index=[
            #    cls.genomic_coordinates({i["ensembl_assembly_name"]: "{0}:{1}:{2}".format(i["ensembl_seq_region_name"],                    i["ensembl_start"],                    i["ensembl_end"])}
            #    ) for i in datalist
            #],
        )

        # extract out any transcripts that are not referenced to current assembly!
        # liftoverneeded  = transcripts.drop(transcripts["ensembl_assembly_name"].str.contains("")

        logger.info("Ensembl: Complete")

        if parameters.store_local:
            cls.save_local(
                transcripts,
                str(thisgene.ensembl_id) + "_" + parameters.active_assembly +
                "_transcripts")

        return transcripts


class Entrez(api):
    # refseq protein link = gene_protein_refseq
    # refseq trancript link = gene_nuccore_refseqrna
    # clinvar link = gene_clinvar
    # together = https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&&db=clinvar&retmode=json&id=672&linkname=gene_clinvar,gene_protein_refseq

    #  can't retrieve an infinite list of ids in 1 go, limit to this number
    retrieval_limit = 500

    # this store the dates of when the different assembly versions were released to work out the exon positions
    assembly_history = None

    def __init__(self):
        raise Exception("This class can not have instances")

    @classmethod
    def esearch(cls, db, term, blockraise=False):

        result = cls.access(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db={0}&retmode=json&term={1}"
            .format(db, term))
        if "ERROR" in result:
            if blockraise:
                print("Warning: Entrez search responded with {0}".format(
                    result["ERROR"]))
            else:
                raise Exception(
                    "Warning: Entrez search responded with {0}".format(
                        result["ERROR"]))
        elif (not "esearchresult" in result
              or len(result["esearchresult"]["idlist"]) == 0):
            if blockraise:
                print(
                    "Warning: Entrez search returned no items for {0} failed".
                    format(term))
            else:
                raise Exception(
                    "Warning: Entrez search responded with {0}".format(
                        result["ERROR"]))

        return result["esearchresult"]["idlist"]

    @classmethod
    def search_and_retrieve(cls, db, term, blockraise=False):
        links = cls.esearch(db, term, blockraise)

        # if there are no results, return back
        if len(links) == 0:
            return {}

        return_var = {}

        for i in range(0, len(links), cls.retrieval_limit):
            data = cls.esummary(db, links[i:i + cls.retrieval_limit])
            # we have data from entrez, now extract
            return_var.update({i: data[i] for i in data["uids"]})

        return (
            return_var  # return the list of transcritps --- currently has no exon data
        )

    @classmethod
    def esummary(cls, db, id):
        """ gets the entrez summary for the supplied id(s) from the supplied db, multiple ids can be provided as CSV or list """

        # do some simple checking
        id = cls.check_ids(id)

        result = cls.access(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db={0}&retmode=json&id={1}"
            .format(db, id))
        try:
            return result["result"]
        except:
            print(result)
            raise Exception("")

    @classmethod
    def efetch(cls, db, id, retmode=None, rettype=None):
        """ gets the entrez entire record for the supplied id(s) from the supplied db, multiple ids can be provided as CSV or list """

        # do some simple checking
        id = cls.check_ids(id)

        return cls.access(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={0}&id={1}&retmode={2}&rettype={3}"
            .format(db, id, retmode, rettype),
            isjson=False,
        )

    @classmethod
    def elink(cls, dbfrom, id, linkname):

        id = cls.check_id(id)
        result = cls.access(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom={0}&retmode=json&id={1}&linkname={2}"
            .format(dbfrom, id, linkname))

        # retrieve the links
        if "ERROR" in result:
            raise Exception("Entrez link responded with {0}".format(
                result["ERROR"]))

        return [
            i["links"] for i in result["linksets"][0]["linksetdbs"]
            if i["linkname"] == linkname
        ][0]

    @classmethod
    def get_transcripts(cls, thisgene):
        """ Provide the Entrez gene id to retrieve RefSeq sequences """
        if type(thisgene) is not gene:
            raise Exception("Function expects gene object")

        cls.check_id(thisgene.entrez_id)

        if parameters.store_local:
            try:
                # try and load local copy

                return_var = cls.load_local(
                    str(thisgene.entrez_id) + "_" +
                    parameters.active_assembly + "_transcripts")
                logger.info("Entrez: Transcript dataset loaded from file.")
                return return_var
            except Exception as ex:
                # if that fails retrieve from api
                pass

        if thisgene.location is None:
            thisgene.load_genomic_location()

        logger.info("Entrez: Retrieving gene table")

        # retrieve the gene table for the gene to extract the genomic coordinates for transcripts / exons
        result = cls.efetch("gene", thisgene.entrez_id, "text",
                            "gene_table").split("\n")
        assembly = ""
        assemblies = cls.get_assemblies()
        ctranscript = ""
        start = 0
        end = 0
        ccds = ""
        tdata = {}
        links = []

        # now process the returned data
        for line in result:
            cline = line.strip()
            # continue if line is empty
            if cline == "":
                continue

            testtranscript = re.search("(?<=RNA)([\s]*N[MR]_[\d]{6,})(.)(\d)",
                                       cline)

            if assembly == "" and "Assembly" in cline:
                # assembly hasn't been set yet and assembly found in this line
                for i in assemblies.index:
                    if i.lower() in cline.lower():
                        # assembly found
                        assembly = i
                        break

                if assembly == "":
                    raise Exception("Assembly not found in Entrez gene table")
            elif testtranscript is not None:
                # outlining the next transcript, so cycle the works to find the id

                ctranscript = testtranscript.group(1).strip()

                if ctranscript not in tdata:
                    tdata[ctranscript] = {}

                tdata[ctranscript].update({
                    "ccdsid": ccds,
                    "exons": [],
                    "start": 0,
                    "end": 0,
                    "strand": 0,
                    "canonical": False
                })

                ccds = ""
                curver = int(testtranscript.group(3).strip())
                transcripthistory = {}
                # here we to work out the exons by getting the version of the trancript that is the most recent compared to the major releases of the assembly accession
                # thisgene.location contains this information
                for v in range(curver, 0, -1):
                    searchresults = cls.search_and_retrieve(
                        "nuccore", "{0}.{1}".format(ctranscript, v))
                    if not len(searchresults) == 1:
                        logger.warn(
                            "Entrez: Expected a single record to be returned by eSearch for {0}.{1}"
                            .format(ctranscript, v))
                        raise Exception(
                            "Expected a single search result for {0}.{1}".
                            format(ctranscript, v))
                    vupdate = datetime.strptime(
                        list(searchresults.values())[0]["updatedate"],
                        "%Y/%m/%d")
                    begindate = thisgene.location[
                        parameters.active_assembly.lower()]["assemblyupdate"]

                    transcripthistory.update({
                        v: {
                            "from": begindate,
                            "to": vupdate,
                            "uid": list(searchresults.keys())[0]
                        }
                    })
                    if v < curver:
                        transcripthistory[v + 1]["from"] = vupdate
                    if vupdate <= thisgene.location[parameters.active_assembly.
                                                    lower()]["assemblyupdate"]:
                        break

                usever = [
                    v["uid"] for k, v in transcripthistory.items()
                    if v["from"] <= thisgene.location[
                        parameters.active_assembly.lower()]["assemblyupdate"]
                    or not (v["from"] > thisgene.location[
                        parameters.active_assembly.lower()]["assemblyupdate"]
                            and v["to"] > thisgene.location[
                                parameters.active_assembly.lower()]
                            ["assemblyupdate"] and not thisgene.location[
                                parameters.active_assembly.lower()]["recent"])
                ][0]

                # add to the links to retrieve the rest of the information later
                links.append(usever)

            elif False and assembly.lower() == parameters.active_assembly:
                # only use this information if the assembly in the gene table is the same as the active assembly.

                # this tabulated data could be used but not yet
                tabulated = [
                    v for v in cline.split("\t") if not v.strip() == ""
                ]

                if len(tabulated) == 1:
                    # still not actual table

                    continue  #
                elif "-" not in tabulated[0]:
                    headers = tabulated

                    continue

                gc = tabulated[0].split("-")
                tdata[ctranscript]["ccdsid"] = ccds

                if not len(headers) == len(tabulated):
                    # table is misformed
                    hascoding = True in [
                        "coding" in v.lower() for v in headers
                    ]

                    if hascoding:
                        # if the transcript has coding headers
                        # this is assumes the misalignment is that this exon has no coding for this transcript
                        # in this case the data column needs to be shifted past the coding column
                        newtabs = []
                        tabindex = 0

                        # count how many "intervals"...should match the "-", if it does this has coding ranges
                        intcount = ["interval" in v.lower()
                                    for v in headers].count(True)
                        dashcount = ["-" in v.lower()
                                     for v in tabulated].count(True)
                        for i, v in enumerate(headers):
                            if "coding" in v.lower() and dashcount < intcount:
                                newtabs.append("")
                            elif i >= len(tabulated):
                                # headers length is longer than the tabulated form
                                newtabs.append("")
                            else:

                                newtabs.append(tabulated[tabindex])
                                tabindex += 1
                        tabulated = newtabs
                    else:
                        # may not have any introns so add on the end
                        for i in range(len(headers) - len(tabulated)):
                            tabulated.append("")

                if gc[0] > gc[1]:
                    # if the start is higher than the end it indicates that it is on the reverse strang
                    tdata[ctranscript]["exons"].append({
                        headers[i].replace(" ", "_").strip(): tabulated[i]
                        for i in range(len(headers))
                    })
                    tdata[ctranscript]["strand"] = -1
                else:
                    tdata[ctranscript]["exons"].append({
                        headers[i].replace(" ", "_").strip(): tabulated[i]
                        for i in range(len(headers))
                    })
                    tdata[ctranscript]["strand"] = 1
                #    tdata[ctranscript]["end"] = gc[1]
                #    if tdata[ctranscript]["start"] == 0:
                #        # only update start if it has no value
                #        tdata[ctranscript]["start"] = gc[0]

                pass

        # the above gene table contains exon information for the latest assembly only!
        # more work is required.
        # first get the transcript link from the gene
        #links = cls.elink("gene", thisgene.entrez_id,"gene_nuccore_refseqrna")

        # if there are no results, return back
        if len(links) == 0:
            return {}

        logger.info("Entrez: Retrieving transcript summaries")
        # retrieve the transcript summaries from entrez
        tdata = {}  # collect everything in dictionary for ease

        for i in range(0, len(links), cls.retrieval_limit):
            data = cls.esummary("nuccore", links[i:i + cls.retrieval_limit])

            # we have data from entrez, now extract
            tdata.update({data[i]["caption"]: data[i]
                          for i in data["uids"]
                          })  # add the data to the tdata collection

        # determine canoncial transcripts.
        # there are 2 ways to do this: either retrieve the genbank file for each transcript which means return of extra data...or query but with potentially invalid result
        logger.info("Entrez: Retrieving transcript records")
        result = cls.efetch("nuccore", list(tdata.keys()), "text", "gb")
        logger.info("Entrez: Parsing data")
        accessions = result.split(
            "//")  # split the returned flat file into individual records
        for acc in accessions:
            if not "KEYWORDS" in acc:
                # empty record
                continue
            accession = re.search("N[MR]_[\d]*",
                                  acc)  # the accession in on the first line
            accession = accession[0] if accession is not None else None
            if accession is None or accession not in tdata:
                continue
            accver = re.search("N[MR]_[\d]*\.\d", acc)  # get the version
            accver = accver[0] if accver is not None else ""
            ccds = re.search("CCDS[\d]*.\d", acc)
            ccds = ccds[0] if ccds is not None else ""
            tdata[accession]["ccdsid"] = ccds
            #exons = re.findall("(?<=exon)[\s]*[\d]*\.\.[\d]*",acc)

            if "RefSeq Select" in acc or "MANE Select" in acc:
                # found canoncial

                tdata[accession]["canonical"] = True

            #if len(exons) > 0:
        #     # has exon data therefore calculate
        #    for exonrange in exons:
        #        eparts = re.findall("[\d]*",exonrange)
        #        # 0 is start 1 is end

        listofdata = [
            cls.rename_keys(v, prefix="entrez_") for v in tdata.values()
        ]

        # entrez start and end positions depend on the strang
        logger.info("Entrez: Complete")

        return_var = pd.DataFrame(data=listofdata)

        if parameters.store_local:
            cls.save_local(
                return_var,
                str(thisgene.entrez_id) + "_" + parameters.active_assembly +
                "_transcripts")

        return return_var
        # return tdata

    @classmethod
    def get_variants(cls, thisgene, search={}, batch=None):
        """ Retrieves clinvar variants for the supplied id """
        if type(thisgene) is not gene:
            raise Exception("Function expects gene object")

        cls.check_id(thisgene.entrez_id)

        if type(batch) is not int:
            batch = cls.retrieval_limit

        # this code attempts to retrieve a local copy of the dataset
        return_var = None

        if parameters.store_local:
            try:
                # try and load local copy

                return_var = cls.load_local(
                    str(thisgene.entrez_id) + "_variants")
                logger.info("Entrez: Variant dataset loaded from file.")
                return return_var
            except Exception as ex:
                # if that fails retrieve from api
                print(ex)
                pass

        if batch > 500:
            logger.info(
                "Entrez doesn't support batches of more than 500 in json mode."
            )
            batch = 500

            # first get the link from the gene

        logger.info("Entrez: Loading from api")

        logger.info("Entrez: Retrieving gene to variant links")

        links = cls.elink("gene", thisgene.entrez_id, "gene_clinvar")

        logger.info("Entrez: {0} links identified".format(len(links)))
        # if there are no results, return back 1
        if len(links) == 0:
            return pd.DataFrame()

        return_var = pd.DataFrame()
        printedone = False
        logger.info(
            "Entrez: Retrieving summaries in batches of {0} ".format(batch))

        for i in range(0, len(links), batch):

            data = cls.esummary("clinvar", links[i:i + batch])

            pagedf = pd.DataFrame([
                cls.rename_keys(data[k], prefix="entrez_")
                for k in data["uids"]
            ])

            # for each of the ids returned by entrez, translate the json into a form pandas can understand and then append to the dataframe
            return_var = pd.concat([return_var, pagedf], ignore_index=True)
            logger.info("Entrez: Progress {0:.2f}%".format(
                min((i + 500) / len(links), 1) * 100))

        # extract the change
        # this is either in the entrez_variantion_set or can use regex to extract from title

        #return_var["entrez_change"] = return_var.apply(lambda x: re.search("c.[^\s]*",x["entrez_title"])[0] if ":c." in x["entrez_title"] else "" ,axis=1)

        return_var = cls.extract_hgvsc(return_var,
                                       "entrez_title",
                                       prefix="entrez_")

        if parameters.store_local:
            cls.save_local(return_var, str(thisgene.entrez_id) + "_variants")

        logger.info("Entrez: Variant dataset load complete.")
        #return_var = return_var.rename(
        #    index=return_var["entrez_title"].apply(cls.variant_coordinates))
        return return_var

    @classmethod
    def check_id(cls, id):
        if (type(id) is str and id == "") or (
            (type(id) is int or str(id).isnumeric())
                and int(id) <= 0) or not (type(id) is str or type(id) is int):
            raise Exception("Function requires valid Entrez ID")
        return id

    @classmethod
    def check_ids(cls, id):
        """ Checks the supplied ID(s). If more than one are supplied, a CSV string will be returned """
        if type(id) is list:
            id = ",".join([cls.check_id(i) for i in id if not i.strip() == ""])
        elif type(id) is int:
            id = cls.check_id(id)
        elif "," in id:
            id = ",".join([
                cls.check_id(i) for i in id.split(",") if not i.strip() == ""
            ])
        return id


class LOVD(api):

    __name = "LOVD"

    # base url https://databases.lovd.nl/shared/api/rest.php
    # /variants/<GENE>?&format=application/json

    #  can't retrieve an infinite list of ids in 1 go, limit to this number
    retrieval_limit = 500

    def __init__(self):
        raise Exception("This class can not have instances")

    @classmethod
    def get_variants(cls, thisgene, raw=False):
        """ If no path is provided, the LOVD api is queried using the recognised gene name (eg BRCA1) and returns a list of associated variants. """
        if type(thisgene) is not gene:
            raise Exception("Function expects gene object")

        if thisgene.symbol == "":
            raise Exception("Valid gene symbol required.")

        # this code attempts to retrieve a local copy of the dataset
        return_var = None
        if parameters.store_local:
            try:
                # try and load local copy

                return_var = cls.load_local(str(thisgene.symbol) + "_variants")
                logger.info("LOVD: Variant dataset loaded from file.")
                return return_var
            except Exception as ex:
                # if that fails retrieve from api
                logger.debug(ex)
                pass

        # if path variable has been set, try load from file...if this fails load from api
        logger.info("LOVD: Retrieving data.")
        result = cls.access(
            "https://databases.lovd.nl/shared/api/rest.php/variants/{0}?format=application/json"
            .format(thisgene.symbol))

        # if there are no results, return back 1
        if len(result) == 0:
            logger.warn("LOVD: No data returned")
            return pd.DataFrame()

        logger.info("LOVD: Parsing data.")
        # serveral parameters in this json are lists but don't have more than one value!
        converted = []

        # first get the active assembly
        active_assembly_ref = cls.get_active_assembly()

        for var in result:

            newvar = {
                cls.rename_keys(k, prefix="lovd_"):
                (cls.parse_data(v[0] if type(v) is list and len(v) == 1 else v)
                 )
                for k, v in var.items()
            }

            converted.append(newvar)

        if not raw:
            # go through the dataset and group by id and variant dna
            return_var = cls.combine_duplicates(
                pd.DataFrame(converted),
                ["lovd_Variant/DBID", "lovd_Variant/DNA"])

            # this may lead to one record having multiple positions
            # the idea is to take the most recently updated.
            return_var = cls.parse_frame(return_var)

            return_var["lovd_grch37_coordinates"] = ""
            return_var["lovd_grch38_coordinates"] = ""

            # extract the variant out

            return_var = cls.extract_hgvsc(return_var,
                                           "lovd_Variant/DNA",
                                           prefix="lovd_")

            return_var["lovd_variant"] = return_var.apply(lambda x: re.sub(
                "NM_[\d]*.\d:c.[\d|+|_|\-|(|)]*", "", x["lovd_Variant/DNA"]),
                                                          axis=1)

            # there maybe issues with the variant in the record...mark issues

            # extract the genomic coordinates
            for index, record in return_var.iterrows():

                positions = record["lovd_position_genomic"]
                updates = record["lovd_edited_date"]

                if type(positions) is not list:
                    positions = [positions]
                if type(updates) is not list:
                    updates = [updates]
                updates = [
                    datetime.strptime(x.replace(":", ""), '%Y-%m-%dT%H%M%S%z')
                    for x in updates
                ]
                setdict = {}

                for i, p in enumerate(positions):
                    for k, v in p.items():
                        this_asssembly = cls.assembly_index(k)[0]
                        update = False
                        conflict = False
                        if "lovd_conflicting_{0}_coordinates".format(
                                this_asssembly) not in list(
                                    return_var.columns):
                            # this assembly hasn't been seen before, add new column
                            return_var["lovd_conflicting_{0}_coordinates".
                                       format(this_asssembly)] = False
                            return_var["lovd_{0}_coordinates".format(
                                this_asssembly)] = ""

                        if k not in setdict:
                            update = True
                        elif not v == setdict[k]["position"]:
                            conflict = True
                            if updates[i] > setdict[k]["lastedit"]:
                                update = True
                            # most recent update

                        if "?" in v or v.strip() == "":
                            conflict = True

                        if conflict:
                            return_var.at[index,
                                          "lovd_conflicting_{0}_coordinates".
                                          format(this_asssembly)] = True

                        if update:
                            setdict[k] = {
                                "position": v,
                                "lastedit": updates[i]
                            }
                            return_var.at[index, "lovd_{0}_coordinates".
                                          format(this_asssembly)] = v
                            update = False

        else:
            return_var = pd.DataFrame(converted)

        if parameters.store_local:
            cls.save_local(return_var, str(thisgene.symbol) + "_variants")

        #return_var = return_var.rename(index=return_var["lovd_Variant/DNA"].
        #                               apply(cls.variant_coordinates))

        # make sure all values are converted out!
        # if the
        #
        logger.info("LOVD: Variant dataset load complete.")
        return return_var


class gene:

    # identifers
    ensembl_id = ""
    entrez_id = ""
    omim_id = ""
    symbol = ""
    ccds_id = ""
    uniprot_id = ""
    hgnc_data = None
    hgnc_id = ""

    location = None

    __transcripts = pd.DataFrame()
    __transcript_assembly = ""
    __variants = pd.DataFrame()
    __variantsources = []

    __api = {}

    def __init__(self,
                 symbol="",
                 entrez_id="",
                 ensembl_id="",
                 omim_id="",
                 api=None):
        """ Initialises a gene. At least one identification is required. Additional APIs can be provided using dictionaries """
        # print(logger.getEffectiveLevel())
        self.ensembl_id = ensembl_id
        self.entrez_id = entrez_id
        self.omim_id = omim_id
        self.symbol = symbol

        self.__api = {
            "hgnc": GeneNames,
            "ensembl": Ensembl,
            "entrez": Entrez,
            "lovd": LOVD,
        }

        if (self.ensembl_id == "" and self.entrez_id == ""
                and self.symbol == "" and self.omim_id == ""):
            raise Exception(
                "Gene initialisation requries at least one identification")
        if type(api) is dict:
            self.__api.update(api)
        elif api is not None:
            raise Exception(
                "Additional APIs must be provided in dictionary form")

        if self.ensembl_id == "" or self.entrez_id == "" or self.symbol == "":
            logger.warn(
                "Some of the identifying information has not been provided, use .get_identification() to retrieve full identification from HGNC API"
            )

        # here get the entrez summary for the gene. this will retrieve the location historys so we can work out the coordinate differences between default_assemblies
        return

    def __str__(self):
        return "Stored information for this gene\nName:\t\t{0} \nEnsembl:\t{1}\nOMIM:\t\t{2}\nEntrez:\t\t{3}\n".format(
            self.symbol, self.ensembl_id, self.omim_id, self.entrez_id)

    def get_identification(self):
        """ Uses the information provided to retrieve the rest from the HGNC API """
        try:
            if not self.symbol == "":

                ids = GeneNames.fetch(self.symbol)["docs"][0]

            elif not self.entrez_id == "":
                ids = GeneNames.fetch(self.entrez_id,
                                      term="entrez_id")["docs"][0]

            elif not self.ensembl_id == "":
                ids = GeneNames.fetch(self.ensembl_id,
                                      term="ensembl_gene_id")["docs"][0]

            elif not self.omim_id == "":
                ids = GeneNames.fetch(self.omim_id, term="omim_id")["docs"][0]

            self.symbol = ids["symbol"]
            self.ensembl_id = ids["ensembl_gene_id"]
            self.entrez_id = ids["entrez_id"]
            self.omim_id = ids["omim_id"]
            self.hgnc = ids
        except Exception as ex:
            raise Exception(
                "Unexpected problem with retrieving identification")

    def load_genomic_location(self):
        """ Accesses the entrez api to retrieve the genomic history across previous assemblies """

        if (type(self.entrez_id) is not int
                and not str(self.entrez_id).isnumeric()) or (
                    type(self.entrez_id) is int and self.entrez_id <= 0):
            logger.warn(
                "Gene has incomplete identification, retrieving information now."
            )
            self.get_identification()

        self.location = {}
        try:
            result = Entrez.esummary("gene",
                                     self.entrez_id)[str(self.entrez_id)]
        except Exception as ex:
            logger.warn(ex)
            raise Exception("The gene with entrez id {0} was not found".format(
                self.entrez_id))
        genomicinfo = result["genomicinfo"][0]
        cur_start = genomicinfo["chrstart"]
        cur_stop = genomicinfo["chrstop"]
        cur_acc = genomicinfo["chraccver"].split(".")[0]
        #cycle through the history of the gene
        #we don't need to use entrez api for each accession verson as the data contains the last major version number in the returned data
        #grch38 is the major version designation. The version in this field represents the first version that is denoted by GRCH38, the previous version will be
        #GRCh37 (the previous major version)
        #This helpful can reduce the number of searches (api calls) we need to do as the location history list is chronological##

        #
        #first limit the list to the chromosome history

        chromohistory = [
            h for h in result["locationhist"]
            if h["chraccver"].split(".")[0] == cur_acc
        ]  #

        skiptover = int(chromohistory[0]["assemblyaccver"].split(".")[1])

        recentassembly = ""
        for i, h in enumerate(chromohistory):
            thisversion = int(chromohistory[i]["assemblyaccver"].split(".")[1])
            #get the current chromohistory data
            if thisversion > skiptover:
                continue  #

            assemblysearch = Entrez.search_and_retrieve(
                "assembly", h["assemblyaccver"])
            if len(assemblysearch) > 1:
                raise Exception(
                    "Search for assembly accession {0} returned multiple results "
                    .format(h["assemblyaccver"]))
            assemblydata = list(assemblysearch.values())[0]

            aname = assemblydata["assemblyname"].split(".")[
                0]  # without the patch number
            assemblies = api.get_assemblies()
            if aname not in list(assemblies.index):
                # add this assembly to the list
                api.add_assembly(aname, [aname])

            self.location.update({
                aname.lower(): {
                    "start":
                    h["chrstart"],
                    "stop":
                    h["chrstop"],
                    "chracc":
                    h["chraccver"],
                    "strand":
                    1 if h["chrstop"] >= h["chrstart"] else -1,
                    "assemblyversion":
                    h["assemblyaccver"],
                    "assemblyupdate":
                    datetime.strptime(assemblydata["asmreleasedate_refseq"],
                                      "%Y/%m/%d %H:%M"),
                    "recent":
                    i == 0,
                }
            })

            #ready the next major release
            skiptover = int(
                assemblydata["lastmajorreleaseaccession"].split(".")[1])
        return self.location

    def load_transcripts(self):
        """ Load the transcripts from Ensembl and Entrez databases, or by from csv that has been previous exported """
        self.__transcripts = pd.DataFrame()

        for k, cur_api in self.__api.items():

            if callable(getattr(cur_api, "get_transcripts", None)):

                # first remove all labels in dataframe that match this key
                # remove all columns with prefix of ensembl_
                self.__transcripts = self.__remove_common_labels(
                    self.__transcripts, "{0}_".format(k.lower()), axis=1)

                # get the data and then merge
                d = cur_api.get_transcripts(self)
                n_ccdsid = [c for c in d.columns if "ccdsid" in c.lower()]
                o_ccdsid = [
                    c for c in self.__transcripts.columns
                    if "ccdsid" in c.lower()
                ]
                # requires a column with ccdsid in the name
                if len(n_ccdsid) == 0:
                    logger.warn(
                        "No column labelled with 'ccdsid' found in dataset from {0}"
                        .format(k))
                elif len(o_ccdsid) == 0:
                    # if no columns are found in the dataset, just overwrite with new
                    self.__transcripts = d
                else:
                    self.__transcripts = self.__transcripts.merge(
                        d,
                        left_on=o_ccdsid[0],
                        right_on=n_ccdsid[0],
                        how="outer")

        if "entrez_uid" in self.__transcripts:
            validator_versions = VariantValidator.get_transcript_versions(
                self.symbol)

            self.__transcripts["validator_version"] = self.__transcripts[
                "entrez_accessionversion"].apply(lambda x: validator_versions[
                    x.split(".")[0]] if type(x) is str and "." in x else "")

        self.__transcript_assembly = parameters.active_assembly
        return self.__transcripts.copy()

    @classmethod
    def __remove_common_labels(cls,
                               df,
                               prefix="",
                               suffix="",
                               anywhere="",
                               axis=0):
        if axis == 0:  # remove index
            thelist = [
                i for i in df.index
                if (i.startswith(prefix) or not prefix == "") and (
                    i.endswith(suffix) or suffix == "") and (
                        anywhere in i or anywhere == "")
            ]
        else:
            thelist = [
                i for i in df.columns
                if (i.startswith(prefix) or prefix == "") and (
                    i.endswith(suffix) or suffix == "") and (
                        anywhere in i or anywhere == "")
            ]

        return df.drop(thelist, axis=axis)

    def transcripts(self,
                    source="entrez",
                    start=-1,
                    end=-1,
                    order=-1,
                    reversed=False):

        if not self.__transcript_assembly == parameters.active_assembly:
            logger.warn(
                "Active assembly is different to transcript loaded data. Reloading."
            )
            self.load_transcripts()

        return self.__transcripts.copy()

    def check_variant_source(self, source=None):

        translate = {"clinvar": "entrez"}

        sources = [
            k.lower() for k, a in self.__api.items()
            if a.get_variants is not None
        ]

        if type(source) is str:
            source = [source]
        if type(source) is not list:
            source = []

        if len(source) == 0:
            source = sources

        source = [
            translate[s.lower()]
            if s.lower() in translate.keys() else s.lower() for s in source
        ]

        # return only valid source
        source = list(set(source).intersection(sources))
        if len(source) == 0:
            raise Exception(
                "Valid source options include a list of one or more from: BX, Entrez, LOVD "
            )
        return source

    def load_variants(self, source=None, path=None):
        """
        Load Variants the variant data set either by api accss or through a file path
        This assumes the data has been directly accessed from the api. 
        Any csv containing mismatch headings will cause issues

        Paths are provided through a dictionary containing the source as the key and path as value.
        Entrez and ClinVar are the same source, so function will default to entrez if both provided
        """

        if path is str:
            # user has provided path, attempt to read and load data
            # TODO check the resulting dataframe to ensure the basic required information is present
            try:
                logging.info("Loading variants from path:")
                self.__variants = api.load_path(path)
                # we don't know what has been loaded so keep this blank
                self.__variantsources = []
                return self.__variants.copy()
            except:
                pass
        elif type(path) is not dict:
            path = {}

        source = self.check_variant_source(source)

        source = source[0]

        for k, cur_api in self.__api.items():
            # if api has correct function
            if callable(
                    getattr(cur_api, "get_variants", None)
            ) and k.lower in source and not k.lower() in self.__variantsources:
                # first remove all labels in dataframe that match this key
                # remove all columns with prefix of ensembl_
                self.__variants = self.__remove_common_labels(self.__variants,
                                                              "{0}_".format(
                                                                  k.lower()),
                                                              axis=1)
                # now add the merge the new data

                #self.__variants = self.__variants.merge(
                #        cur_api.get_variants(self),
                #        left_index=True,
                #        right_index=True,
                #        how="outer")

                self.__variantsources.append(k.lower())
        #  validate the source list

        return self.__variants.copy()


# useful functions
def print_columns(df, step=3, padding=50):
    """ Prints the columns of the supplied dataframe """
    cur = 0
    for c in df.columns:
        columnsneeded = math.floor(len(c) / padding)
        print(c.ljust(padding + (padding * columnsneeded), " "), end="")
        cur += 1 + columnsneeded
        if cur >= step:
            cur = 0
            print("")
