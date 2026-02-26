import sys
import time
import subprocess
import os  
from multiprocessing import Pool, RLock

def main():
    if len(sys.argv) != 3:
        sys.exit(f'python {sys.argv[0]} <script_path> <num_cpus>\n')

    script_path = sys.argv[1]
    num_cpus = int(sys.argv[2])

    log_file_name, datetime_str = get_datetime_strings(script_path)
    init_log_file(log_file_name, datetime_str)

    lock = RLock()
    with Pool(processes=num_cpus, initializer=init, initargs=(lock,)) as pool:
        with open(script_path) as script_file:
            for line in script_file:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                pool.apply_async(run_script_line, (line, log_file_name))

        pool.close()
        pool.join()

    finalize_log_file(log_file_name, datetime_str)

def run_script_line(script_line, log_file_name):
    global lock  # Declare lock as global
    try:
        result = subprocess.run(script_line, shell=True, stderr=subprocess.PIPE, check=True)
        success = True
    except subprocess.CalledProcessError as e:
        result = e
        success = False

    current_time = time.localtime()
    timestamp = f'{current_time.tm_year}-{current_time.tm_mon:02}-{current_time.tm_mday:02} {current_time.tm_hour:02}:{current_time.tm_min:02}:{current_time.tm_sec:02}'
    lock.acquire()
    try:
        with open(log_file_name, 'a') as log_file:  # Use 'a' to append to the log file
            if success:
                log_file.write(f'{timestamp} Success -> {script_line}\n')
            else:
                log_file.write(f'{timestamp} Failure -> {script_line}\n')
                log_file.write(f'Error: {result.stderr.decode()}\n')
    finally:
        lock.release()

def get_datetime_strings(script_path):
    current_time = time.localtime()
    datetime_str = f'{current_time.tm_year}{current_time.tm_mon:02}{current_time.tm_mday:02}{current_time.tm_hour:02}{current_time.tm_min:02}{current_time.tm_sec:02}'
    log_file_name = f'{os.path.splitext(script_path)[0]}.{datetime_str}.plog'
    return log_file_name, f'{current_time.tm_year}-{current_time.tm_mon:02}-{current_time.tm_mday:02} {current_time.tm_hour:02}:{current_time.tm_min:02}:{current_time.tm_sec:02}'

def init(lock_param):
    global lock
    lock = lock_param

def init_log_file(log_file_name, datetime_str):
    with open(log_file_name, 'a') as log_file:  # Use 'a' to append to the log file
        log_file.write(f'=== START {datetime_str} START ===\n')

def finalize_log_file(log_file_name, datetime_str):
    with open(log_file_name, 'a') as log_file:
        log_file.write(f'=== END {datetime_str} END ===\n')

if __name__ == '__main__':
    main()