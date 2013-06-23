//
//  SoapBubbleViewController.m
//  SoapBubble
//
//  Created by Kevin Utecht on 6/18/13.
//

#import "SoapBubbleViewController.h"
#import "RayTracer.h"



@interface SoapBubbleViewController ()
@property (weak, nonatomic) IBOutlet UISlider *thicknessSlider;
@property (weak, nonatomic) IBOutlet UIImageView *mainImageView;
@property (weak, nonatomic) IBOutlet UIImageView *previewImageView;
@property (weak, nonatomic) IBOutlet UIBarButtonItem *renderButton;

@property BOOL mainInProgress;
@property BOOL previewInProgress;
@end

@implementation SoapBubbleViewController


- (void)viewDidLoad
{
    [super viewDidLoad];
    self.previewInProgress = NO;
}

- (void)renderMainView:(UIBarButtonItem*)sender
{
    if (!self.mainInProgress)
    {
        UIActivityIndicatorView *spinner = [[UIActivityIndicatorView alloc] initWithActivityIndicatorStyle:UIActivityIndicatorViewStyleGray];
        [spinner startAnimating];
        self.navigationItem.rightBarButtonItem = [[UIBarButtonItem alloc] initWithCustomView:spinner];
        dispatch_queue_t mainRenderQ = dispatch_queue_create("Render Main Queue", NULL);
        dispatch_async(mainRenderQ, ^{
            UIImage *tracedImage = raytrace(self.mainImageView.bounds.size.width,
                                            self.mainImageView.bounds.size.width,
                                            self.thicknessSlider.value);
            dispatch_async(dispatch_get_main_queue(), ^{
                self.navigationItem.rightBarButtonItem = sender;
                self.mainImageView.image = tracedImage;
            });
        });
    }
}

- (void)renderPreviewView
{
    if (!self.previewInProgress)
    {
        self.previewInProgress = YES;
        dispatch_queue_t previewRenderQ = dispatch_queue_create("Render PreviewQueue", NULL);
        dispatch_async(previewRenderQ, ^{
            UIImage *tracedImage = raytrace(self.previewImageView.bounds.size.width,
                                            self.previewImageView.bounds.size.width,
                                            self.thicknessSlider.value);
            dispatch_async(dispatch_get_main_queue(), ^{
                self.previewImageView.image = tracedImage;
                self.previewInProgress = NO;
            });
        });
    }
}

- (void)viewDidAppear:(BOOL)animated
{
    [super viewDidAppear:animated];
    [self renderPreviewView];
}

- (void)didReceiveMemoryWarning
{
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

- (IBAction)render:(UIBarButtonItem *)sender {
    [self renderMainView:sender];
}

- (IBAction)sliderValueChanged:(id)sender {
    [self renderPreviewView];
}

@end
