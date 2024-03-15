//// Generate cross products that keep track of both variable names and values

def crossProduct(Map<String, List<?>> mapOfLists, boolean dryRun) {
    if (dryRun) {
        def result = [:]
        for (key in mapOfLists.keySet()) {
            def list = mapOfLists[key]
            result[key] = list.first()
        }
        return Channel.of(result)
    } else 
        crossProduct(mapOfLists)
}

def crossProduct(Map<String, List<?>> mapOfLists) {
    def keys = new ArrayList(mapOfLists.keySet())
    def list = _crossProduct(mapOfLists, keys)
    return Channel.fromList(list)
}

def _crossProduct(mapOfLists, keys) {
    if (keys.isEmpty()) {
        return [[:]]
    }
    def key = keys.remove(keys.size() - 1)
    def result = []
    for (recursiveMap : _crossProduct(mapOfLists, keys))
        for (value : mapOfLists.get(key)) {
            def copy = new LinkedHashMap(recursiveMap)
            copy[key] = value 
            result.add(copy)
        }
    return result
}


//// Collect CSVs from many execs

def collectCSVs(inputChannel) { collectCSVsProcess(inputChannel.map{toJson(it)}.toList()) }

def toJson(tuple) {
    def args = tuple[0]
    def path = tuple[1]
    return groovy.json.JsonOutput.toJson([
        path: path.toString(), 
        args: args
    ])
}

process collectCSVsProcess {
    debug true
    cache true
    scratch false
    time { 10.m * Math.pow(2, task.attempt-1) }
    memory { 4.GB * Math.pow(2, task.attempt-1) }
    errorStrategy 'retry'
    maxRetries '4'
    input: 
        // we use file not path, as we want the json strings to be dumped into files
        file jsonObjects
    output:
        path aggregated
    publishDir deliverables(workflow), mode: 'copy', overwrite: true
    """
    aggregate
    """
}



////

process setupPigeons {
    debug true
    executor 'local'
    scratch false
    time { 2.h * Math.pow(2, task.attempt-1) }
    errorStrategy 'retry'
    maxRetries '4'
    input:
        env JULIA_DEPOT_PATH
        path julia_env
    output:
        path julia_env
    script:
        template 'setup_pigeons.sh'
}


/////

process head {
    debug true 
    executor 'local'
    input:
        path aggregated
    output:
        path aggregated
    """
    # TODO: column command is not supported by bash inside container
    for i in `ls $aggregated/*.csv`
    do
        echo \$i
        head \$i | sed 's/,,/,NA,/g' | column -t -s, 
    done
    """
}

//// Git utils

process checkGitUpdated {
    cache false 
    debug true
    executor 'local'
    input:
        val dryRun
    shell:
    if (dryRun) " " else
    '''
    cd !{projectDir}
    git remote update
    UPSTREAM=${1:-'@{u}'}
    LOCAL=$(git rev-parse @)
    REMOTE=$(git rev-parse "$UPSTREAM")
    BASE=$(git merge-base @ "$UPSTREAM")

    if [ $LOCAL = $REMOTE ]; then
        echo "Branch $(git rev-parse --abbrev-ref HEAD) up to date"
        exit 0
    elif [ $LOCAL = $BASE ]; then
        echo "Need to pull"
    elif [ $REMOTE = $BASE ]; then
        echo "Need to push"
    else
        echo "Diverged"
    fi
    git status
    '''
}

process commit {
    executor 'local'
    input:
        path token // pass a channel produced at very end to make sure we commit only at completion 
        val dryRun
    script:
    if (dryRun) " " else
    """
    cd ${projectDir}
    git add deliverables/*
    git commit -m "Auto: add deliverables for run ${workflow.runName} of ${workflow.scriptName}"
    git push
    """
}


/////

def pow(int i, int j) { java.lang.Math.pow(i, j) as Integer}

////

def deliverables(workflow) { 'deliverables/' + workflow.scriptName.replace('.nf','') }
