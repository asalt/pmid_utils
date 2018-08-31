import os
import sqlite3 as sql
import time
import re

import bs4
from bs4 import BeautifulSoup
from Bio import Entrez, Medline
from bcmproteomics import ispec
import numpy as np
import pandas as pd

from functools import lru_cache

email = 'saltzman@bcm.edu'
if email is None:
    print('Put in your email!')
    email = input('> ')

Entrez.email = email

__db__ = 'PMID.db'
__basedir__ = os.path.expanduser('~')

sql.register_adapter(np.int64, lambda val: int(val))
sql.register_adapter(np.int32, lambda val: int(val))

class OpenDB(object):
    """
    Simple CM for sqlite3 databases. Commits everything at exit.
    """
    def __init__(self, path=None):

        if path is None:
            path = os.path.join(__basedir__, __db__)

        self.path = path
        self._conn = None
        self._cursor = None

    @property
    def conn(self):
        if self._conn is None:
            self._conn = sql.connect(self.path, check_same_thread=False)
        return self._conn

    @property
    def cursor(self):
        if self._cursor is None:
            self._cursor = self.conn.cursor()
        return self._cursor


    def __enter__(self):
        self._conn = sql.connect(self.path, check_same_thread=False)
        self._cursor = self.conn.cursor()
        self.create()
        return self

    def __exit__(self, exc_class, exc, traceback):
        self.conn.commit()
        self.conn.close()


    def create(self):

        publication_def = '''
        CREATE TABLE IF NOT EXISTS publication(
        PubMedID INTEGER PRIMARY KEY,
        title TEXT,
        abstract TEXT,
        author TEXT,
        terms TEXT,
        pubdate TEXT
        )
        '''

        gene2pubmed_def = '''
        CREATE TABLE IF NOT EXISTS gene2pubmed(
        ID INTEGER PRIMARY KEY AUTOINCREMENT,
        TaxonID INT NOT NULL,
        GeneID INT NOT NULL,
        PubMedID INT NOT NULL,
        FOREIGN KEY(PubMedID) REFERENCES publication(id)
        UNIQUE(TaxonID, GeneID, PubMedID)
        )
        '''

        pubmed_genecount = '''
        CREATE TABLE IF NOT EXISTS pubmed_genecount(
        PubMedID INTEGER PRIMARY KEY,
        GeneCount INT NOT NULL,
        FOREIGN KEY(PubMedID) REFERENCES publication(id)
        )
        '''

        self.conn.execute(publication_def)
        self.conn.execute(gene2pubmed_def)
        self.conn.execute(pubmed_genecount)
        self.conn.commit()


    # def add(self, geneid, pubmed_id, taxonid):
    def add(self, taxonid_gid_pmid_collection, download_abstracts=False):
        '''
        [ taxonid, geneid, pmid,
          taxonid, geneid, pmid,
          ...
        ]
        '''

        # for row in taxonid_gid_pmid_collection:
        #     cursor.execute("""
        #     INSERT INTO gene2pubmed(TaxonID, GeneID, PubMedID)
        #     """)

        self.conn.executemany(
            """INSERT OR IGNORE INTO gene2pubmed(TaxonID, GeneID, PubMedID)
            VALUES (?, ?, ?)""",
            taxonid_gid_pmid_collection
        )
        self.conn.commit()


        if download_abstracts:
            pubmed_ids = [x[2] for x in taxonid_gid_pmid_collection]
            self.fetch_publications(pubmed_ids)

    def get_publications(self, pubmed_ids):
        """
        get publications
        """
        params = ', '.join(('?')*len(pubmed_ids))
        df = pd.read_sql("""
        SELECT * from publication where PubMedID in ({})
        """.format(params), con=self.conn, params=(*pubmed_ids,))
        return df

    def fetch_publications(self, pubmed_ids, geneid):
        """
        fetch publications related to a geneid
        """
        cursor = self.conn.execute('select PubMedID from publication where PubMedID in (%s)' %
                                ','.join('?'*len(pubmed_ids)), pubmed_ids)
        existing_pmids = [x for y in cursor.fetchall() for x in y]

        missing_pmids = [x for x in pubmed_ids if int(x) not in existing_pmids]

        if missing_pmids:
            pmid_info = self.fetch_pmid_info(missing_pmids)
            self.update_publications(pmid_info)

        params = ', '.join(('?')*len(pubmed_ids))
        df = pd.read_sql("""select gene2pubmed.TaxonID, gene2pubmed.GeneID, gene2pubmed.PubMedID,
        publication.title, publication.abstract, publication.author, publication.terms, publication.pubdate
        from gene2pubmed
        INNER JOIN publication on gene2pubmed.PubMedID = publication.PubMedID
        where publication.PubMedID IN ({}) AND gene2pubmed.GeneID = ?
        """.format(params),
                         con=self.conn, params=(*pubmed_ids, geneid), index_col='PubMedID')

        return df



    def fetch_pmid_info(self, pmids):
        pmids = [str(x) if not isinstance(x, str) else x for x in pmids]

        print('Querying NCBI for abstracts')
        with Entrez.efetch(db='pubmed', id=pmids, rettype='medline') as handle:
            # medline is easy to parse!
            records = list(Medline.parse(handle))  # returns a list of dictionaries
            # time.sleep(3)
            time.sleep(.3)
        pmid_info = {record.get('PMID'): record for record in records}
        time.sleep(.3)
        # self.count_genes(pmid_info)
        return pmid_info

    @staticmethod
    def get_collection(value, key='AU'):
        """utility function for returning either a joined list (as str) or empty string"""
        ret = value.get(key)
        if ret and isinstance(ret, (list, tuple)):
            return '; '.join(ret)
        elif isinstance(ret, str):
            return ret
        else:
            return ''

    def update_publications(self, pmid_info):

        'title abstract author terms'
        print('Saving {} abstracts'.format(len(pmid_info)))

        insertions = [
            (key, value.get('TI', ''), value.get('AB', ''), self.get_collection(value, 'AU'),
             self.get_collection(value, 'OT'), value.get("DP", ''))
            for key, value in pmid_info.items()
        ]

        self.conn.executemany(
            """INSERT OR IGNORE INTO publication(PubMedID, title, abstract, author, terms, pubdate)
            VALUES(?, ?, ?, ?, ?, ?)
            """, insertions
        )
        self.conn.commit()

    def get_pmids(self, geneid, thresh=50):

        cursor = self.conn.execute("""SELECT gene2pubmed.PubMedID from gene2pubmed
        INNER JOIN pubmed_genecount on gene2pubmed.PubMedID = pubmed_genecount.PubMedID
        where GeneID = ? and pubmed_genecount.GeneCount < ?
        """, (geneid, thresh))

        return [x for y in cursor.fetchall() for x in y]

    def load_gene2pmid(self):
        #TODO: download from NCBI and update
        #TODO: more than just human

        import csv
        gene2pubmed = os.path.join(os.path.split(__file__)[0], 'gene2pubmed')

        data = list()
        with open(gene2pubmed) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row['#tax_id'] != '9606':
                    continue
                data.append((row['#tax_id'], row['GeneID'], row['PubMed_ID']))

        self.add(data)

    def count_gids_and_update(self):
        self.cursor.execute

        df = pd.read_sql("""
        SELECT PubMedID, count(PubMedID) from gene2pubmed
        GROUP BY PubMedID
        """, self.conn).rename(columns={'count(PubMedID)': 'GeneCount'})
        df.to_sql('pubmed_genecount', self.conn, if_exists='replace', index=False)


    def count_total(self, geneid, thresh=50):

        cursor = self.cursor.execute("""
        SELECT COUNT(DISTINCT gene2pubmed.PubMedID)
        from gene2pubmed
        inner join pubmed_genecount on gene2pubmed.PubMedID = pubmed_genecount.PubMedID
        where GeneID = ? and pubmed_genecount.GeneCount < ?
        """, (geneid, thresh)
        )

        return cursor.fetchall()[0][0]

    @lru_cache()
    def get_and_fetch(self, geneid, thresh):
        pmids = self.get_pmids(geneid, thresh=thresh)
        df = self.fetch_publications(pmids, geneid)
        return df


    def keyword_filter(self, geneid, keywords, case=False, pmid_gene_cutoff=50):

        # if isinstance(case, int):
        #     if not case:
        #         flags = repeat(re.I, len(keywords))
        #     else:
        #         flags = repeat(0, len(keywords))
        # else:
        #     flags = [re.I if not bool(c) else 0 for c in case]

        # kws = [re.compile(kw, flag) for kw, flag in zip(keywords, flags)]

        # pmids = self.get_pmids(geneid, thresh=pmid_gene_cutoff)
        # df = self.fetch_publications(pmids, geneid)

        df = self.get_and_fetch(geneid, pmid_gene_cutoff)

        boolean = list()
        regx = '|'.join(keywords)
        for col in ('title', 'abstract', 'author', 'terms'):
            df[col] = df[col].fillna('')
            bool_res = df[col].str.contains(regx, case=case, regex=True)
            boolean.append(bool_res)
        boolean_frame = pd.concat(boolean, axis=1).any(1)

        return df[boolean_frame]