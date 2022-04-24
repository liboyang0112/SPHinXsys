#!/bin/bash
if [ -z ${SPH_DIR+x} ] ; then
	SOURCE="${BASH_SOURCE[0]}"
	while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
	  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
	  SOURCE="$(readlink "$SOURCE")"
	  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
	done
	export SPH_DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
	alias sphcd='cd $SPH_DIR'
fi