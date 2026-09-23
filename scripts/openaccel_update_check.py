#!/usr/bin/env python3
"""Validate ordered public update captures and pack actual reference outputs."""
import argparse
import hashlib
import io
import json
from pathlib import Path
from openaccel_reference_check import CONTRACT, integer, numbers, require, unique_object

WIDTHS = (3,1,10,24,9,25,31,10,6,2)
OUTPUTS = (1,1,3,1,1,1,6,2,1,1)
PHASES = (2,2,5,3,3,3,4,1,1,1)


def close(a,b):
    return abs(a-b) <= 2e-12*max(1,abs(a),abs(b))


def validate(stage,x,y):
    def flags(values):
        require(all(v in (0,1) for v in values),'invalid flag')
    if stage==0:
        require(0 < x[2] <= 1,'invalid pressure relaxation')
    elif stage==2:
        flags(x[9:10])
    elif stage in (3,5):
        require(x[0]>0 and 0 < x[23] <= 1,'invalid flux scales')
        if stage==5: flags(x[24:25])
    elif stage==4:
        require(x[0]>0 and 0 < x[8] <= 1,'invalid inlet scales')
    elif stage==6:
        flags(x[3:6]+x[30:31]+y[3:6])
        require(x[3]==x[4]==x[5] and y[3]==y[4]==y[5],'partial face reversal flag')
        require(sum(sum(x[21+3*s+j] for s in range(3))**2 for j in range(3))>0,'zero face normal')
    elif stage==7:
        flags(x[9:10])
        if x[9]==0:
            require(close(sum(x[3:6]),1) and sum(a*a for a in x[6:9])>0,'bad trace quadrature')
    elif stage==8:
        flags(x[5:6]); require(0<=x[3]<=1,'invalid trace blend')
    elif stage==9:
        require(x[1]>0,'no open outlet area: pinned trace mean is undefined')


def load_updates(directory):
    records,hashes={},{}
    iterations=set()
    for path in sorted(directory.glob('*.jsonl')):
        content=path.read_bytes()
        rows=[json.loads(line,object_pairs_hook=unique_object) for line in content.splitlines()]
        require(len(rows)>=10,'empty/truncated update file')
        h,t=rows[0],rows[-1]
        require(h.get('kind')=='header' and h.get('schema')==1 and h.get('fixture')=='public_channel'
                and h.get('producer')=='openaccel' and h.get('ranks')==1,'not a single-rank public update capture')
        require(h.get('reference_revision')==CONTRACT['reference']['revision']
                and h.get('solver_revision')==CONTRACT['reference']['solver_gitlink']['revision'],'wrong source pin')
        iteration=integer(h['iteration'],'iteration')
        require(iteration not in iterations,'duplicate iteration'); iterations.add(iteration)
        phase=-1; count=0; trace_stage=0
        for row in rows[1:-1]:
            if row.get('kind')=='phase':
                require(type(row.get('value')) is int and row['value']==phase+1 and row['value']<=7,'wrong update order')
                phase=row['value']; continue
            require(row.get('kind')=='update' and row.get('phase')==phase,'invalid update event')
            stage=integer(row['stage'],'stage',0)
            require(stage<10 and PHASES[stage]==phase,'update uses wrong SIMPLE phase')
            if stage in (7,9,8):
                position={7:0,9:1,8:2}[stage]
                require(position>=trace_stage and position<=trace_stage+1,'wrong trace update order')
                trace_stage=position
            entity=integer(row['entity'],'entity'); sample=integer(row['sample'],'sample',0)
            require(sample < (6 if stage==3 else 3 if stage in (4,5,7,8) else 1),'invalid sample')
            key=(stage,iteration,entity,sample); require(key not in records,'duplicate update identity')
            if stage==1:
                require((0,iteration,entity,0) in records,'pressure increment stored before pressure update')
            x,y=row['inputs'],row['outputs'];numbers(x,WIDTHS[stage],'update inputs');numbers(y,OUTPUTS[stage],'update outputs')
            validate(stage,x,y); records[key]=row; count+=1
        require(phase==7 and t.get('kind')=='end' and type(t.get('records')) is int and t['records']==count,'incomplete update sequence')
        hashes[path.name]=hashlib.sha256(content).hexdigest()
    require(iterations and iterations==set(range(1,max(iterations)+1)),'missing iteration')
    for it in iterations:
        groups=[{(n,s):r for (a,i,n,s),r in records.items() if (a,i)==(stage,it)} for stage in range(10)]
        require(all(groups),'missing update stage')
        require(groups[0].keys()==groups[1].keys()==groups[2].keys(),'pressure/increment/velocity node coverage differs')
        require(groups[5].keys()==groups[7].keys()==groups[8].keys(),'outlet flux/trace coverage differs')
        for stage,samples in ((3,6),(4,3),(5,3)):
            ids={n for n,s in groups[stage]}
            require(set(groups[stage])=={(n,s) for n in ids for s in range(samples)},'missing entity sample')
        require(set(groups[6])=={(n,0) for n,s in groups[5]},'outlet reversal faces differ')
        require(set(groups[9])=={(1,0)},'one public outlet patch required')
        mean=groups[9][(1,0)]
        require(close(sum(r['outputs'][0] for r in groups[7].values()),mean['inputs'][0])
                and close(sum(r['outputs'][1] for r in groups[7].values()),mean['inputs'][1]),'trace mean uses different moments')
        for (n,s),r in groups[5].items():
            reversal=groups[6][(n,0)]
            require(r['outputs'][0]==reversal['inputs'][s],'reversal does not consume stored flux')
            require(r['inputs'][24]==reversal['inputs'][3+s],'flux used a different reversal flag')
            trace=groups[8][(n,s)]
            require(trace['inputs'][2]==mean['outputs'][0],'trace uses a different patch mean')
            require(groups[7][(n,s)]['inputs'][9]==trace['inputs'][5]==r['inputs'][24],'trace and flux flags differ')
            require(trace['outputs'][0]==reversal['inputs'][18+s],'reversal uses a different trace')
            if it>1:
                previous=records[(6,it-1,n,0)]
                require(r['inputs'][22]==previous['outputs'][s],'outlet flux history was changed')
                require(r['inputs'][24]==previous['outputs'][3+s],'outlet flag history was changed')
        if it>1:
            for stage,offset in ((3,22),(4,7)):
                for (n,s),r in groups[stage].items():
                    require(r['inputs'][offset]==records[(stage,it-1,n,s)]['outputs'][0],'stored mass flux history differs')
    return records,hashes


def pack_updates(directory,output):
    records,hashes=load_updates(directory)
    require(not output.exists() and not output.with_suffix('.json').exists(),'output/provenance exists')
    stream=io.StringIO();stream.write('MARS_PUBLIC_UPDATE_REPLAY_V1 {}\n'.format(len(records)))
    for key,row in sorted(records.items()):
        stream.write(' '.join(str(v) for v in key)+'\n')
        stream.write(' '.join(format(v,'.17g') for v in row['inputs'])+'\n')
        values=row['outputs']+[0]*(6-len(row['outputs']))
        stream.write(' '.join(format(v,'.17g') for v in values)+'\n')
    with output.open('x') as handle: handle.write(stream.getvalue())
    output.with_suffix('.json').write_text(json.dumps(dict(records=len(records),source_sha256=hashes,
        packed_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),coverage='ordered-frozen-updates'),indent=2)+'\n')
    print('Prepared {} public update records: {}'.format(len(records),output))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('exports',type=Path);parser.add_argument('--pack',type=Path)
    args=parser.parse_args()
    try:
        if args.pack: pack_updates(args.exports,args.pack)
        else:
            records,hashes=load_updates(args.exports)
            print('PASS: update capture order/coverage/history: {} records, {} files; numerical replay pending'.format(len(records),len(hashes)))
    except (ValueError,KeyError,TypeError,OSError,OverflowError) as error:
        parser.exit(1,'ERROR: '+str(error)+'\n')
