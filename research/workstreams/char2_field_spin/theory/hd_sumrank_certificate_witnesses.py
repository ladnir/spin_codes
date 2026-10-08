"""Pinned rational witnesses for the s=16 feedback-refresh proof candidate.

Only verified all-word intervals belong here. The final replay independently
checks both coverage and the summed failure bound; incomplete lists are not
advertised as a certificate. These are proof inputs, not performance data.
"""

PAIRED=(
    (1,1,'992087/993431','1','656/486783',8),
    (2,2,'212815/213392','1','697/258275',8),
    (3,3,'976001/979973','1','4000/988167',8),
    (4,4,'466491/469024','1','4731/876872',8),
    (5,5,'781952/787263','1','5299/786098',8),
    (6,7,'253191/255256','1','3125/386534',8),
    (8,9,'295428/298645','1','5343/496250',8),
    (10,11,'663875/672388','1','12157/960579',8),
    (12,14,'823589/837078','1/2','1053/129698',8),
    (15,17,'911788/930493','1/2','6685/659338',8),
    (18,22,'797093/817862','1/4','4175/610762',8),
    (23,28,'120061/124026','1/4','7603/881234',8),
    (29,34,'802000/835243','1/6','7979/996627',8),
    (35,41,'617911/648772','1/8','339/41867',12),
    (42,50,'575083/610380','1/8','1600/161781',12),
    (51,60,'930051/999236','1/10','4282/401403',12),
    (61,72,'690383/752868','2/25','3707/315779',12),
    (73,84,'851189/943431','2/25','5342/383445',12),
    (85,97,'521493/588269','2/25','9267/568408',12),
    (98,110,'298697/342927','7/100','10444/586543',12),
    (111,123,'724931/848201','7/100','13381/661209',12),
    (124,136,'363501/433450','7/100','671/29605',12),
    (137,149,'803307/974896','3/50','5296/223113',12),
    (150,163,'637236/788149','3/50','22979/882214',12),
    (164,180,'223857/283676','3/50','28014/967723',12),
    (181,198,'348963/453203','3/50','14383/450472',12),
    (199,217,'420169/560049','3/50','23455/668631',16),
    (218,237,'724755/992908','3/50','18051/470104',16),
    (238,258,'317966/448375','3/50','33673/803974',16),
    (259,281,'412607/599746','3/50','44118/968941',16),
    (282,306,'366771/551125','3/50','33936/685345',16),
    (307,336,'294657/459040','1/20','36221/708196',16),
    (337,371,'609019/990784','1/20','49073/874977',16),
    (372,413,'553724/947527','1/20','53456/864293',16),
    (414,464,'225181/409421','1/20','38467/559407',16),
    (465,527,'56375/110334','1/20','59925/776783',16),
    (528,604,'391635/839488','1/20','4895/55973',20),
    (605,699,'409350/980659','1/20','96731/967568',20),
    (700,816,'327447/899794','1/20','98758/856213',20),
    (817,962,'16263/52915','1/20','62038/462521',20),
)

GLOBAL=(
    (1024,'138671/581931','39308/782933','39366/242489'),
    (1220,'91846/506861','33543/663773','113732/591217'),
    (1400,'83721/594728','10071/197351','193659/880562'),
    (1536,'100168/861357','19988/388321','161021/669049'),
    (1780,'42223/511658','48066/917309','236055/850049'),
    (2048,'23891/421437','45311/843271','239700/753677'),
    (2500,'10572/352643','54352/956227','96955/251008'),
    (3072,'9937/731904','24688/398175','375986/799811'),
    (4096,'1672/554579','55742/688735','274418/437001'),
    (4500,'1310/800951','70213/776569','266531/386259'),
    (6144,'95/709988','135341/777887','885955/966397'),
    (7400,'1/65533','140333/140340','458505/458512'),
)


def replay(progress=None):
    from flint import arb,ctx
    from hd_sumrank_paired_certificate import receipt as paired
    from hd_sumrank_global_certificate import receipt as global_bound
    from hd_randomized_outer_proof import dense_tail_receipt
    with ctx.workprec(256):
        covered=set(); total=arb(0); receipts=[]
        for first,last,z,alpha,p,cutoff in PAIRED:
            item=paired(first,last,z,alpha,p,cutoff=cutoff)
            item.pop('rejected_rows')
            receipts.append(item)
            if progress: progress(item)
        for center,z,alpha,p in GLOBAL:
            item=global_bound(center,z,alpha,p)
            if not item['accepted']: raise ArithmeticError('Pinned global witness did not certify')
            receipts.append(item)
            if progress: progress(item)
        for item in receipts:
            for first,last in item['certified_intervals']: covered.update(range(first,last+1))
            if item['union_log2_upper'] is not None:
                total+=arb(2)**arb(item['union_log2_upper'])
        tail=dense_tail_receipt()
        covered.update(range(tail['active_rows_first'],tail['active_rows_last']+1))
        total+=arb(2)**arb(tail['union_log2_upper'])
        missing=[r for r in range(1,8193) if r not in covered]
        bound=total.log()/arb(2).log()
        return dict(construction='sumrank-s16-independent-feedback-refresh',message_length=2**20,
                    code_length=2**21,field_size=65536,rate='1/2',forbidden_maximum_weight=838860,
                    state=16,packet_width=64,uncovered_active_rows=missing,union_log2_upper=str(bound),
                    directed=True,whole_code_certificate=not missing and bool(bound < -40),
                    native_implementation=False,performance_claim=False,receipts=len(receipts)+1)


if __name__=='__main__':
    import json
    print(json.dumps(replay(lambda row:print(json.dumps(row),flush=True))),flush=True)
