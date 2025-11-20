import os
import shutil
from datetime import datetime

def backup_file(filepath):
    """Create a backup of a file with timestamp."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
        
    # Create backups directory if it doesn't exist
    backup_dir = os.path.join(os.path.dirname(filepath), 'backups')
    os.makedirs(backup_dir, exist_ok=True)
    
    # Create backup filename with timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    filename = os.path.basename(filepath)
    backup_name = f"{os.path.splitext(filename)[0]}_{timestamp}{os.path.splitext(filename)[1]}"
    backup_path = os.path.join(backup_dir, backup_name)
    
    # Copy file to backup
    shutil.copy2(filepath, backup_path)
    print(f"Created backup: {backup_path}")
    
def restore_backup(backup_path, target_path=None):
    """Restore a file from backup."""
    if not os.path.exists(backup_path):
        raise FileNotFoundError(f"Backup not found: {backup_path}")
        
    # If no target specified, restore to original location
    if target_path is None:
        target_path = os.path.join(
            os.path.dirname(os.path.dirname(backup_path)),
            os.path.basename(backup_path).split('_')[0] + os.path.splitext(backup_path)[1]
        )
    
    # Create backup of current file before restoring
    if os.path.exists(target_path):
        backup_file(target_path)
        
    # Restore backup
    shutil.copy2(backup_path, target_path)
    print(f"Restored from backup: {backup_path} to {target_path}")
    
def list_backups(filepath):
    """List all backups for a given file."""
    backup_dir = os.path.join(os.path.dirname(filepath), 'backups')
    if not os.path.exists(backup_dir):
        print("No backups found")
        return []
        
    filename = os.path.basename(filepath)
    base_name = os.path.splitext(filename)[0]
    backups = [f for f in os.listdir(backup_dir) if f.startswith(base_name)]
    backups.sort(reverse=True)  # Most recent first
    
    if backups:
        print("\nAvailable backups:")
        for backup in backups:
            print(f"- {backup}")
    else:
        print("No backups found")
        
    return backups
