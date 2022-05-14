#! /usr/bin/env python
import screed
import sourmash
from sourmash.index import Index


class ScreedIndex(Index):
    def __init__(self, filename, selection_dict=None):
        self.filename = filename
        if selection_dict:
            self.selection_dict = selection_dict
        else:
            self.selection_dict = dict(abund=False)

    def __len__(self):
        return 1

    @property
    def location(self):
        return self.filename

    def signatures(self):
        sd = self.selection_dict
        assert 'ksize' in sd
        assert 'moltype' in sd
        assert 'scaled' in sd or 'num' in sd
        assert 'abund' in sd

        if sd.get('containment'):
            assert sd.get('num')

        ksize = sd['ksize']
        moltype = sd['moltype']
        scaled = sd.get('scaled', 0)
        num = sd.get('num', 0)
        track_abundance = sd.get('abund', False)

        is_protein = False
        dayhoff = False
        hp = False
        if moltype == 'DNA':
            pass
        elif moltype == 'protein':
            is_protein=True
        elif moltype == 'dayhoff':
            dayhoff = True
        elif moltype == 'hp':
            hp = True

        mh = sourmash.MinHash(num, ksize, is_protein=is_protein,
                              dayhoff=dayhoff, hp=hp,
                              track_abundance=track_abundance,
                              scaled=scaled)
                              
        for record in screed.open(self.filename):
            this_mh = mh.copy_and_clear()
            this_mh.add_sequence(record.sequence, force=True)

            ss = sourmash.SourmashSignature(this_mh, name=record.name)
            yield ss

    def select(self, **kwargs):
        selection_dict = dict(self.selection_dict)
        for k, v in kwargs.items():
            if k in selection_dict:
                if selection_dict[k] != v:
                    raise ValueError(f"cannot select on two different values for {k}")
            selection_dict[k] = v

        return ScreedIndex(self.filename, selection_dict)

    def insert(self, signature):
        raise NotImplemented

    def save(self, *args, **kwargs):
        raise NotImplemented

    def load(self, *args, **kwargs):
        raise NotImplemented
