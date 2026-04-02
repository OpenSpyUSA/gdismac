#ifdef __APPLE__
#import <Cocoa/Cocoa.h>
#include "mac_activation.h"

void mac_activate_foreground_app_now(void)
{
@autoreleasepool
  {
  NSApplication *app = [NSApplication sharedApplication];

  if ([app respondsToSelector:@selector(setActivationPolicy:)])
    [app setActivationPolicy:NSApplicationActivationPolicyRegular];

  [[NSRunningApplication currentApplication]
      activateWithOptions:NSApplicationActivateAllWindows];
  [app activateIgnoringOtherApps:YES];
  }
}

gboolean mac_activate_foreground_app(gpointer data)
{
(void) data;

mac_activate_foreground_app_now();

return(FALSE);
}
#endif
